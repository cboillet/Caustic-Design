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

private:
    //std::vector<Point> Xo;
    std::vector<double> c_weights;
    std::vector<Vertex_handle> p_vertices;
    //std::map<Point, Vertex_handle> map;


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


    void runInterpolation();
    void runInterpolation(QString imageFile, QString datFile);
    std::vector<Point>& getXo(){return Xo;}
    void computeWeights(std::vector<Vertex_handle> neighbors, Point oP);
    void findNaturalNeighbor(Point oP);
    void findNeighbor(Point oP);
    //std::vector<Point> computeXr(std::vector<std::pair<Vertex_handle, FT> > vertices_weight);



};

void multiThreadInterpolation(Interpolation* inter, uint id, uint from, uint to);
void findNeighborr(Point oP, Interpolation* inter, uint id);

#endif // INTERPOLATION_H

