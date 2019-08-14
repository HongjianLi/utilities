#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

void write(const vector<string>& lines, const string filename)
{
	ofstream ofs(filename);
	for (const auto& line : lines)
	{
		ofs << line << endl;
	}
}

int main(int argc, char* argv[])
{
	const string delimiter = "@<TRIPOS>MOLECULE";
	vector<string> lines;
	size_t id = 0;
	for (string line; getline(cin, line);)
	{
		if (line == delimiter && !lines.empty())
		{
			write(lines, to_string(++id) + ".mol2");
			lines.clear();
		}
		lines.push_back(line);
	}
	write(lines, to_string(++id) + ".mol2");
//	cout << "Splitted into " << id << " files" << endl;
}
