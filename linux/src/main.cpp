#include "h/CMesh.h"
#include <iostream>
#include <cstring>
#include <ctime>
#include <cstdlib>
using namespace std;

#define ON 1
#define OFF 0
#define DEBUG ON

int main(int arg, char *args[]) {
    srand((unsigned)time(NULL));
    if (arg < 2) {
        cout << "Error: Need Input Arguments!" << endl;
        return -1;
    }
    string filename(args[1]);
    size_t found = filename.rfind(".obj");
    if (found == string::npos) {
        cout << "Error: File Format Is Not [.obj]!" << endl;
        return -1;
    }
    string writefile = filename.substr(0, found) + ".iDT.obj";
    CMesh *mesh = new CMesh();
    mesh->Load(args[1]);
    mesh->iDT();
    mesh->Save(writefile.c_str());
    delete mesh;
    return 0;
}
