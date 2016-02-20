#include <iostream>
#include <fstream>
#include <math.h>

using namespace::std;
int main(int argc, char* argv[])
{
	double a = 7.0;
	double b = 0.1;

	double x1 = 0;
	double x2 = 0;
	double x3 = 0;
	double f;

	ifstream fin;
	fin.open("input.dat");
	fin >> x1 >> x2 >> x3;
	fin.close();

	f = sin(x1) + a*pow(sin(x2), 2.0) + b * pow(x3, 4.0) * sin(x1);
	
	ofstream fout;
	fout.open("output.dat");
	fout << f << endl;
	fout.close();
}
