#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

int main(int argc, char* argv[])
{
	const string delimiter = "@<TRIPOS>MOLECULE";
	string id, line;
	size_t numConformers;
	for (bool delimited = false; getline(cin, line);) {
		if (delimited) {
			delimited = false;
			if (line == id) {
				++numConformers;
			} else {
				if (id.size()) {
					if (numConformers == 1) {
						cout << id << endl;
					} else {
						cout << id << '\t' << numConformers << endl;
					}
				}
				id = line;
				numConformers = 1;
			}
		} else if (line == delimiter) {
			delimited = true;
		}
	}
	if (id.size()) {
		if (numConformers == 1) {
			cout << id << endl;
		} else {
			cout << id << '\t' << numConformers << endl;
		}
	}
}
