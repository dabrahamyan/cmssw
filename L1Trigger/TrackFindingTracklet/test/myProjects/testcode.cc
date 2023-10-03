#include <iostream>
#include <string>
#include <fstream>

using namespace std;

int main () {
    bool a = "false";

    a = a << 7;

    cout << a << endl;

    return 0;
}