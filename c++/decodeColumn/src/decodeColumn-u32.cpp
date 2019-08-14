#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint32_t g; // typedef unsigned int uint32_t
	while (cin.read(reinterpret_cast<char*>(&g), sizeof(g)))
	{
		cout << g << endl;
	}
}
