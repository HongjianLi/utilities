#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	float g;
	const auto addr = reinterpret_cast<char*>(&g);
	const auto size = sizeof(g);
	for (string line; getline(cin, line);)
	{
		g = stof(line);
		cout.write(addr, size);
	}
}
