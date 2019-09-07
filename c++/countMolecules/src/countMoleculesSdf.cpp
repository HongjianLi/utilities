#include <iostream>
#include <string>
using namespace std;

int main(int argc, char* argv[])
{
	const string delimiter = "$$$$";
	string id, line;
	size_t numConformers;
	for (bool delimited = true; getline(cin, line);) {
		if (delimited) {
			delimited = false;
			if (line == id) {
				++numConformers;
			} else {
				if (id.size()) {
					cout << id << '\t' << numConformers << endl;
				}
				id = line;
				numConformers = 1;
			}
		} else if (line == delimiter) {
			delimited = true;
		}
	}
	if (id.size()) {
		cout << id << '\t' << numConformers << endl;
	}
}
