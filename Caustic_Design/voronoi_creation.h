#ifndef VORONOI_CREATION_H
#define VORONOI_CREATION_H

#include <scene.h>
#include <QtCore>
#include <QString>

class Scene;
class VoronoiCreator
{
public:
    VoronoiCreator(Scene* sc){}
    ~VoronoiCreator(){}

    bool generate_voronoi(Scene* sc, unsigned npoints, double epsilon);

    void init_points(int npoints, Scene* sc);
    void apply_lloyd_optimization(Scene* sc);

};

#endif // VORONOI_CREATION_H
