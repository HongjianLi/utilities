#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint16_t g; // typedef unsigned short int uint16_t
	const auto addr = reinterpret_cast<char*>(&g);
	const auto size = sizeof(g);
	while (cin.read(addr, size))
	{
		cout << g << endl;
	}
}
