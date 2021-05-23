#include <iostream>
#include <fstream>
using namespace std;
int main(void)
{
    double eyeX, eyeY, eyeZ;
    double lookX, lookY, lookZ;
    double upX, upY, upZ;
    double fovY, aspectRatio, near, far;
    ifstream scene;
    scene.open("scene.txt");
    scene >> eyeX >> eyeY >> eyeZ;
    scene >> lookX >> lookY >> lookZ;
    scene >> upX >> upY >> upZ;
    scene >> fovY >> aspectRatio >> near >> far;

    string cmd;
    while(true)
    {
        scene >> cmd;
        if (cmd == "triangle")
        {
            cout << "triangle";
        }
        else if (cmd == "translate")
        {
            double tx, ty, tz;
            scene >> tx >> ty >> tz;
            cout << "translate";
        }
        else if (cmd == "scale")
        {
            double sx, sy, sz;
            scene >> sx >> sy >> sz;
            cout << "scale";
        }
        else if (cmd == "rotate")
        {
            cout << "rotate";
        }
        else if (cmd == "push")
        {
            cout << "push";
        }
        else if (cmd == "pop")
        {
            cout << "pop";
        }
        else if (cmd == "end")
        {
            break;
        }
    }

    scene.close();

    //cout << eyeX << eyeY<<eyeZ <<"  " <<lookX<<lookY<<lookZ<<endl;
    //cout << upX<<upY<<upZ << "  " << fovY<<aspectRatio<<near<<far << endl;
    return 0;
}
