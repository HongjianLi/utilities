#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint8_t g; // typedef unsigned char uint8_t
	for (string line; getline(cin, line);)
	{
		g = stoul(line);
		cout.write(reinterpret_cast<char*>(&g), sizeof(g));
	}
}
