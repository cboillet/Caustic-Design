#ifndef INTERPOLATION_H
#define INTERPOLATION_H
// STL
#include <map>
#include <vector>

//Qt
#include <QString>

// local
#include "matrix/sparse_matrix.h"
#include "types.h"
#include "scene.h"
#include "window.h"

class Scene;
class VoronoiCreator;
class MainWindow;

class Interpolation{

private:
    std::vector<Point> Xo;
    std::vector<double> c_weights;
    std::vector<Vertex_handle> p_vertices;

public:
    Scene* source_scene;
    Scene* m_scene;
    Scene* compute_scene;
    MainWindow* win;

    Interpolation(Scene* sc, Scene* tsc, Scene* csc, MainWindow* win);
    ~Interpolation(){}

    void runInterpolation();

    std::vector<Point>& getXo(){return Xo;}
    void computeWeights(std::vector<Vertex_handle> neighbors, Point oP);
    void findNaturalNeighbor(Point oP);
    //std::vector<Point> computeXr(std::vector<std::pair<Vertex_handle, FT> > vertices_weight);



};

#endif // INTERPOLATION_H

