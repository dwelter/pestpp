#include <iostream>
#include "RunManagerExternal.h"

using namespace std;

RunManagerExternal::RunManagerExternal(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, const string &_ext_filename, int _max_n_failure)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
	_insfile_vec, _outfile_vec, stor_filename, _max_n_failure), ext_filename(_ext_filename)
{
	cout << "              starting external run manager ..." << endl << endl;
}

void RunManagerExternal::run()
{
	vector<int> waiting_run_ids = get_outstanding_run_ids();

	//if there are outstanding runs that need to be made
	//exit so the external run manager can be involked
	if (!waiting_run_ids.empty())
	{
		ofstream fout_ext(ext_filename);
		fout_ext << get_run_filename() << endl;
		fout_ext << max_n_failure << endl;
		fout_ext.close();

		cout << endl << endl;
		cout << "    External run mananager involked.  Leaving PEST++ ..." << endl;
		cout.flush();
		exit(0);
	}

}

RunManagerExternal::~RunManagerExternal()
{
}
