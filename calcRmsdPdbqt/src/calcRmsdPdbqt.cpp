#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
using namespace std;

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cout << "calcRmsdPdbqt reference.pdbqt docked.pdbqt" << endl;
		return 1;
	}

	string line;
	vector<array<double, 3>> ref;
	for (ifstream ifs(argv[1]); getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			const string element = line.substr(77, 2);
			if (!(element == "H " || element == "HD"))
			{
				ref.push_back
				({
					stod(line.substr(30, 8)),
					stod(line.substr(38, 8)),
					stod(line.substr(46, 8)),
				});
			}
		}
	}
	const double n_inv = 1.0 / ref.size();

	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(6);
	double se = 0;
	size_t i = 0;
	for (ifstream ifs(argv[2]); getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			const string element = line.substr(77, 2);
			if (!(element == "H " || element == "HD"))
			{
				const double dx = ref[i][0] - stod(line.substr(30, 8));
				const double dy = ref[i][1] - stod(line.substr(38, 8));
				const double dz = ref[i][2] - stod(line.substr(46, 8));
				se += dx * dx + dy * dy + dz * dz;
				++i;
			}
		}
		else if (record == "TORSDO")
		{
			cout << sqrt(se * n_inv) << endl;
			se = i = 0;
		}
	}
}
