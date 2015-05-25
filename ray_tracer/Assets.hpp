#ifndef ASSETS_H
#define ASSETS_H

#include <string>

using namespace std;


namespace Assets{
    extern string modelPath;

    extern string texturePath;

    extern string models[];
    extern string textures[];

    extern int modelAmount;
    extern int textureAmount;

    void update();
}


#endif // ASSETS_H
