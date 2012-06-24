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

#include <iostream>
#include <fstream>
#include "MorrisMethod.h"
#include <Eigen/Dense>


using namespace std;
using Eigen::MatrixXd;



int main(int argc, char* argv[])
{
	cout << "Starting Program" << endl;

	MatrixXd b_star_mat;

	MorrisMethod morris(8);

	b_star_mat = morris.create_P_star_mat(7);
	cout << b_star_mat << endl << endl;

	cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
	char buf[256];
    gets(buf);
}