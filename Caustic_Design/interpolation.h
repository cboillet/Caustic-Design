#ifndef INTERPOLATION_H&
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

class Scene;
<<<<<<< HEAD
class VoronoiCreator;

=======
>>>>>>> 13c14360646e9cad8f42118e91dc7f290592b60b
class Interpolation{

private:
    std::vector<Point> Xo;
    std::vector<double> c_weights;
    std::vector<Point> Xr;
    std::vector<Vertex_handle> p_vertices;

public:
<<<<<<< HEAD
    Scene* m_scene;
    Scene* target_scene;
    Scene* compute_scene;

    Interpolation(Scene* sc, Scene* tsc, Scene* csc);
    ~Interpolation(){}

    void runInterpolation();
    bool prepareData();

    std::vector<Point>& getXo(){return Xo;}
    std::vector<Vertex_handle> findNaturalNeighbor(Point oP);
    std::vector<std::pair<Point, FT> > computeWeights(std::vector<Point> neighbors, Point oP){}
=======
    Interpolation(Scene* sc);
    ~Interpolation(){}

    std::vector<Point>& getXo(){return Xo;}

    std::vector<Point> findNaturalNeighbor(Point oP,Scene* sc);
    std::vector<std::pair<Point, FT> > computeWeights(std::vector<Point> neighbors, Point oP,Scene* sc){}
>>>>>>> 13c14360646e9cad8f42118e91dc7f290592b60b
    std::vector<Point> computeXr(std::vector<std::pair<Point, FT> > cWeights, std::vector<Vertex_handle> p_vertices,Scene* sc){}

    std::vector<Vertex_handle> compareCell(Vertex_handle vc);


};

#endif // INTERPOLATION_H

