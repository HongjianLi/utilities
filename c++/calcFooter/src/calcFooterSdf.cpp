#include <iostream>
using namespace std;
int main()
{
	for (string line; getline(cin, line);)
	{
		if (line == "$$$$")
		{
			cout << cin.tellg() << endl;
		}
	}
}
