#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint8_t g; // typedef unsigned char uint8_t
	const auto addr = reinterpret_cast<char*>(&g);
	const auto size = sizeof(g);
	while (cin.read(addr, size))
	{
		cout << g << endl;
	}
}
