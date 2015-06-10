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
    std::vector<Point> Xr;
    std::vector<Vertex_handle> p_vertices;

public:
    Scene* m_scene;
    Scene* target_scene;
    Scene* compute_scene;
    MainWindow* win;

    Interpolation(Scene* sc, Scene* tsc, Scene* csc, MainWindow* win);
    ~Interpolation(){}

    void runInterpolation();

    std::vector<Point>& getXo(){return Xo;}
    std::vector<Vertex_handle> findNaturalNeighbor(Point oP);
    std::vector<std::pair<Vertex_handle, FT> > computeWeights(std::vector<Vertex_handle> neighbors, Point oP);
    std::vector<Point> computeXr(std::vector<std::pair<Vertex_handle, FT> > vertices_weight);



};

#endif // INTERPOLATION_H

