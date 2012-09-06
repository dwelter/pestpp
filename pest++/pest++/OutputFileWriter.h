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
#ifndef OUTPUTFILEWRITER_H
#define OUTPUTFILEWRITER_H

#include <string>
#include <iostream>
class Observations;
class ObjectiveFunc;
class Parameters;
class TranOffset;
class TranScale;
class Jacobian;
class QSqrtMatrix;
class ParameterGroupInfo;

class OutputFileWriter
{
public:
	static void write_rei(const std::string &filename, int iter_no, const Observations &obs,
		const Observations &sim, const ObjectiveFunc &obj_func, const Parameters &pars);
	static void write_par(const std::string &filename, const Parameters &pars, const TranOffset &offset_tran, const TranScale &scale_tran);
	static void write_sen_header(std::ostream &fout, const std::string &case_name);
	static void append_sen(std::ostream &fout, int iter_no, const Jacobian &jac, const ObjectiveFunc &obj_func, const ParameterGroupInfo &par_grp_info);
};
#endif /* OUTPUTFILEWRITER_H */