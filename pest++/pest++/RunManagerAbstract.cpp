/*  
    © Copyright 2012, David Welter
    
    This file is part of PEST++.
   
    PEST++ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PEST++ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/
#include "RunManagerAbstract.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>
#include <cstring>
#include "Transformable.h"
#include "ModelRunPP.h"
#include "Pest.h"
#include "utilities.h"



RunManagerAbstract::RunManagerAbstract(const ModelExecInfo &_model_exec_info)
	 : nruns(0), total_runs(0)
{
	comline_vec = _model_exec_info.comline_vec;
	tplfile_vec = _model_exec_info.tplfile_vec;
	inpfile_vec =_model_exec_info.inpfile_vec;
	insfile_vec = _model_exec_info.insfile_vec;
	outfile_vec =_model_exec_info.outfile_vec;
}


 Observations RunManagerAbstract::get_obs_template(double value) const
 {
	Observations ret_obs;

	for(int i=0, nobs=obs_name_vec.size(); i<nobs; ++i)
	{
		ret_obs[obs_name_vec[i]] = value;
	}
	return ret_obs;
 }