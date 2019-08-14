#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint32_t g; // typedef unsigned int uint32_t
	const auto addr = reinterpret_cast<char*>(&g);
	const auto size = sizeof(g);
	for (string line; getline(cin, line);)
	{
		g = stoul(line);
		cout.write(addr, size);
	}
}
