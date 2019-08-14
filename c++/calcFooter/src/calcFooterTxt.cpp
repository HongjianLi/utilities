#include <iostream>
using namespace std;
int main()
{
	for (string line; getline(cin, line);)
	{
		cout << cin.tellg() << endl;
	}
}
