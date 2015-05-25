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

class Interpolation{

private:
    std::vector<Point> Xo;
    std::vector<double> c_weights;
    std::vector<Point> Xr;
    std::vector<Vertex_handle> p_vertices;
    Scene* sc;

public:
    Interpolation();
    ~Interpolation(){delete sc;}

    Scene getSceneCp();
    Scene& getScene();

    std::vector<Point> findNaturalNeighbor(Point oP);
    std::vector<std::pair<Point, FT> > computeWeights(std::vector<Point> neighbors, Point oP){}
    std::vector<Point> computeXr(std::vector<std::pair<Point, FT> > cWeights, std::vector<Vertex_handle> p_vertices){}



};

#endif // INTERPOLATION_H

