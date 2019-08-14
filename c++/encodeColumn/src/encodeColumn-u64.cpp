#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint64_t g; // typedef unsigned long int uint64_t
	const auto addr = reinterpret_cast<char*>(&g);
	const auto size = sizeof(g);
	for (string line; getline(cin, line);)
	{
		g = stoul(line);
		cout.write(addr, size);
	}
}
