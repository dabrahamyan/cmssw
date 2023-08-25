#include <iostream>
#include <string>
#include <fstream>

using namespace std;

int main () {
    int myArray[3][3];
    
    int counter = 0;
    for (int i = 0; i < 3; i++) {
        for (int ii = 0; ii < 3; i++) {
            counter++;
            myArray[i][ii] = counter;
        }
    }

    cout << sizeof(myArray)/4 << endl;

    return 0;
}