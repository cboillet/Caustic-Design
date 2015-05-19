#ifdef WITH_BOOST
    #include <boost/filesystem.hpp>
#endif

#include "Assets.hpp"

#include <string>

using namespace std;

int Assets::modelAmount = 21;
string Assets::models[] = {"bialetti.off",
                   "eight.off",
                   "killeroo.off",
                   //"plane4x4.off",
                           "quad.off",
                   "bunny.off",
                   "europemap.off",
                   "sphere.off",
                   "bunnysimple.off",
                   "dice.off",
                   "hand.off",
                   "lightcycle.off",
                   "teapot.off",
                   "camel_head.off",
                   "dragon.off",
                   "heptoroid.off",
                   "lucy.off",
                   "torus_tri.off",
                   "cow.off",
                   "drei.off",
                   "homer.off",
                   "mannequin.off"};

int Assets::textureAmount = 19;
string Assets::textures[] = {"checker2.ppm",
                     "checker.ppm",
                     "flowers.ppm",
                     "waterbump.ppm",
                     "world.ppm",
                     "checker3.ppm",
                     "earthbump1k_2.ppm",
                     "grid_normals.ppm",
                     "wood.ppm",
                     "worldspec.png",
                     "checker4.ppm",
                     "earthbump1k.ppm",
                     "grid.ppm",
                     "worldbump.ppm",
                     "worldspec.ppm",
                     "checker5.ppm",
                     "fine.ppm",
                     "lines.ppm",
                     "world.jpg"};

string Assets::modelPath = "meshes/";
string Assets::texturePath = "data/";

void Assets::update(){
    #ifdef WITH_BOOST

    #else
    #endif
}
