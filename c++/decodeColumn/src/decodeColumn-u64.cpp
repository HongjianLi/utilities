#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint64_t g; // typedef unsigned long int uint64_t
	while (cin.read(reinterpret_cast<char*>(&g), sizeof(g)))
	{
		cout << g << endl;
	}
}
