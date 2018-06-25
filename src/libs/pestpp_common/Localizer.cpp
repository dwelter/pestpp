#include <random>
#include <iomanip>
#include <unordered_set>
#include <iterator>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "ParamTransformSeq.h"
#include "ObjectiveFunc.h"
#include "RedSVD-h.h"
#include "covariance.h"
#include "PerformanceLog.h"
#include "system_variables.h"
#include "Localizer.h"

bool Localizer::initialize(PerformanceLog *performance_log)
{
	stringstream ss;
	string filename = pest_scenario_ptr->get_pestpp_options().get_ies_localizer();
	if (filename.size() == 0)
		return false;
	string ext = filename.substr(filename.size() - 3, 3);
	pest_utils::upper_ip(ext);
	if ((ext == "JCB") || (ext == "JCO"))
	{
		performance_log->log_event("loading localizer from binary file");
		mat.from_binary(filename);
	}
	else if (ext == "MAT")
	{
		performance_log->log_event("loading localizer from binary file");
		mat.from_ascii(filename);
	}
	else
	{
		ss << "unrecognnized localizer extension '" << ext << "', should be JCB, JCO, or MAT";
		performance_log->log_event("error: "+ss.str());
		throw runtime_error(ss.str());
	}


	//error checking and building up container of names
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
	set<string> par_names(names.begin(), names.end());
	names = pest_scenario_ptr->get_ctl_ordered_obs_names();
	set<string> obs_names(names.begin(), names.end());

	map<string, vector<string>> pargp_map;
	ParameterGroupInfo *pi = pest_scenario_ptr->get_base_group_info_ptr();

	for (auto &pg : pest_scenario_ptr->get_ctl_ordered_par_group_names())
	{
		
		names.clear();
		for (auto &p : par_names)
			if (pi->get_group_name(p) == pg)
				names.push_back(p);
		pargp_map[pg] = names;
	}

	map<string, vector<string>> obgnme_map;
	for (auto &og : pest_scenario_ptr->get_ctl_ordered_obs_group_names())
	{
		ObservationInfo *oi = pest_scenario_ptr->get_observation_info_ptr();
		names.clear();
		for (auto &o : obs_names)
			if (oi->get_group(o) == og)
				names.push_back(o);
		obgnme_map[og] = names;
	}

	vector<string> missing, dups;
	vector<vector<string>> obs_map;
	set<string> dup_check;
	for (auto &o : mat.get_row_names())
	{
		if (obs_names.find(o) != obs_names.end())
		{
			obs_map.push_back(vector<string>{o});
			if (dup_check.find(o) != dup_check.end())
				dups.push_back(o);
			dup_check.emplace(o);
		}
		else if (obgnme_map.find(o) != obgnme_map.end())
		{
			obs_map.push_back(obgnme_map[o]);
			for (auto &oo : obgnme_map[o])
			{
				if (dup_check.find(oo) != dup_check.end())
					dups.push_back(oo);
				dup_check.emplace(oo);

			}
		}
		else
			missing.push_back(o);
	}
	if (missing.size() > 0)
	{
		ss << " the following rows in " << filename << " were not found in the observation names or observation group names: ";
		for (auto &m : missing)
			ss << m << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}
	if (dups.size() > 0)
	{
		ss << " the following observations were listed more than once (possibly through an obs group): ";
		for (auto & d : dups)
			ss << d << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}

	vector<vector<string>> par_map;
	dup_check.clear();
	for (auto &p : mat.get_col_names())
	{
		if (par_names.find(p) != par_names.end())
		{
			par_map.push_back(vector<string>{p});
			if (dup_check.find(p) != dup_check.end())
				dups.push_back(p);
			dup_check.emplace(p);
		}
		else if (pargp_map.find(p) != pargp_map.end())
		{
			par_map.push_back(pargp_map[p]);
			for (auto &pp : pargp_map[p])
			{
				if (dup_check.find(pp) != dup_check.end())
					dups.push_back(pp);
				dup_check.emplace(pp);
			}
		}
		else
			missing.push_back(p);
	}
	if (missing.size() > 0)
	{
		ss << " the following cols in " << filename << " were not found in the active parameter names or parameter group names: ";
		for (auto &m : missing)
			ss << m << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}
	if (dups.size() > 0)
	{
		ss << " the following parameters were listed more than once (possibly through a par group): ";
		for (auto & d : dups)
			ss << d << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}


	//map all the nz locations in the matrix
	map<int, vector<int>> idx_map;
	int i, j;
	for (int k = 0; k < mat.e_ptr()->outerSize(); ++k)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(*mat.e_ptr(), k); it; ++it)
		{
			if (idx_map.find(it.row()) == idx_map.end())
				idx_map[it.row()] = vector<int>{ it.col() };
			else
				idx_map[it.row()].push_back(it.col());
			//std::cout << "(" << it.row() << ","; // row index
			//std::cout << it.col() << ")\t"; // col index (here it is equal to k)

		}
	}

	//populate the localizer map
	vector<string> vobs, vpar;
	
	for (auto &idx : idx_map)
	{
		vobs = obs_map[idx.first];
		vpar.clear();
		for (auto &i : idx.second)
		{
			vpar.insert(vpar.end(), par_map[i].begin(), par_map[i].end());
		}
		pair<vector<string>, vector<string>> p(vobs,vpar);
		localizer_map.push_back(p);
	}

	cout << "done" << endl;
}