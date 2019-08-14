#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint32_t g; // typedef unsigned int uint32_t
	for (string line; getline(cin, line);)
	{
		g = stoul(line);
		cout.write(reinterpret_cast<char*>(&g), sizeof(g));
	}
}
