#ifndef INTERPOLATION_H
#define INTERPOLATION_H
// STL
#include <map>
#include <unordered_map>
#include <vector>
#include <thread>
#include <mutex>

//Qt
#include <QString>

// local
#include "matrix/sparse_matrix.h"
#include "types.h"
#include "scene.h"
#include "window.h"
#include "glviewer.h"
#include "console_color.h"
#include "timer.h"

class Scene;
class VoronoiCreator;
class MainWindow;

class Interpolation{

public:
    Scene** computeScenes;
    Scene* source_scene;
    Scene* m_scene;
    Scene* compute_scene;
    MainWindow* win;
    std::vector<Point> Xo;
    std::map<Point, Vertex_handle> map;

    Interpolation(Scene* sc, Scene* tsc, Scene* csc, int sitesAmount, MainWindow* win);
    ~Interpolation(){}


    void runInterpolation(QString imageFile, QString datFile);

};

void multiThreadInterpolation(Interpolation* inter, uint id, uint from, uint to);
void findNeighborr(Point oP, Interpolation* inter, uint id, int index);

#endif // INTERPOLATION_H

