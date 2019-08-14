#include <iostream>
using namespace std;
int main()
{
	for (string line; getline(cin, line);)
	{
		if (line[0] == 'T')
		{
			cout << cin.tellg() << endl;
		}
	}
}
