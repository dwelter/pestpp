// This example problem was taken from Saltelli, etel, "Sensitiyivity Analysis in Practice", 2000
// It is present on pg88 and is reproduced from Morris's 1991 paper.

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{
	double b = 0;
	double b_i[20];
	double b_ij[20][20];
	double x[20];
	double w[20];
	ofstream fout;
	ifstream fin;

	fout.open("morris_1991.out");
	fin.open("morris_1991.inp");
	int i;
	for (i = 0; i < 20 && fin; ++i)
	{
		fin >> x[i];
	}
	fin.close();
	if (i < 20) exit(0);


	std::default_random_engine generator;
	generator.seed(1);
	std::normal_distribution<double> distribution(0.0, 1.0);

	


	for (int i = 0; i < 20; ++i)
	{
		if (i == 2 || i == 4 || i == 6)
		{
			w[i] = 2.0 * (1.1 * x[i] /(x[i] +.1) - .5);
		}
		else
		{
			w[i] = 2.0 * (x[i] - .5);
		}
	}


	double y = 0;
	for (int i = 0; i < 20; ++i)
	{
		b_i[i] = (i < 10) ? 20 : distribution(generator);
		for (int j = 0; j < 20; ++j)
		{
			if (j < 6)
			{
				b_ij[i][j] = -15.0;
			}
			else
			{
				b_ij[i][j] = distribution(generator);
			}
		}
	}

	y = 0;
	for (int i = 0; i < 20; ++i)
	{
		y += b_i[i] * w[i];
		for (int j = i+1; j < 20; ++j)
		{
			y += b_ij[i][j] * w[i]*w[j];
			for (int l = j+1; l < 5; ++l)
			{
				y += -10.0 * w[i] * w[j] * w[l];
				for (int s = l+1; s < 4; ++s)
				{
					y += 5.0 * w[i] * w[j] * w[l] * w[s];
				}
			}
		}
	}
	fout << y << endl;
	cout << "Morris_ex complete ..." << endl;
	fout.close();

}