#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	float g;
	const auto addr = reinterpret_cast<char*>(&g);
	const auto size = sizeof(g);
	while (cin.read(addr, size))
	{
		cout << g << endl;
	}
}
