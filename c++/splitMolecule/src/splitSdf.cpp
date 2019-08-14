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
	const string delimiter = "$$$$";
	vector<string> lines;
	size_t id = 0;
	for (string line; getline(cin, line);)
	{
		lines.push_back(line);
		if (line == delimiter)
		{
			write(lines, to_string(++id) + ".sdf");
			lines.clear();
		}
	}
	if (lines.size() > 2)
	{
		write(lines, to_string(++id) + ".sdf");
	}
	cout << "Splitted into " << id << " files" << endl;
}
