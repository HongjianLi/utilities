#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	float g;
	for (string line; getline(cin, line);)
	{
		g = stof(line);
		cout.write(reinterpret_cast<char*>(&g), sizeof(g));
	}
}
