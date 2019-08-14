#include <iostream>
using namespace std;

int main(int argc, char* argv[])
{
	uint8_t g; // typedef unsigned char uint8_t
	while (cin.read(reinterpret_cast<char*>(&g), sizeof(g)))
	{
		cout << g << endl;
	}
}
