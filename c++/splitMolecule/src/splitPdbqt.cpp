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
	vector<string> lines;
	bool content = false;
	size_t id = 0;
	for (string line; getline(cin, line);)
	{
		if (line[0] == 'R')
		{
			content = true;
		}
		if (content)
		{
			lines.push_back(line);
		}
		if (line[0] == 'T')
		{
			content = false;
			write(lines, to_string(++id) + ".pdbqt");
			lines.clear();
		}
	}
	cout << "Splitted into " << id << " files" << endl;
}
