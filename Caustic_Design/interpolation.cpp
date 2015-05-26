#include "interpolation.h"
#include "random.h"

Interpolation::Interpolation(Scene* sc){
//    sc = new Scene();
    if (!sc->getDomain().is_valid()) return;
    double dx = sc->getDomain().get_dx();
    double dy = sc->getDomain().get_dy();
    while (Xo.size() != 500)
    {
        double x = random_double(-dx, dx);
        double y = random_double(-dy, dy);
        Xo.push_back( Point(x, y) );
    }
}

//Scene Interpolation::getSceneCp(){return *sc;}
//Scene& Interpolation::getScene(){return (*sc);}

/*
Method 1:
    - insert point
    - recompute voronoi Diagram
    - see which vertices have changed
 to improve: high computation time!

Method 2:
    Use Delaunay triangulation http://cgl.uni-jena.de/pub/Workshops/WebHome/cgl12new.pdf
    Compute neighbirs using CGAL librairy

Method 3:
    Kirkpatrick method http://cgm.cs.mcgill.ca/~athens/cs507/Projects/2002/PaulSandulescu/
    O(logn)

Current: method 1
*/

std::vector<Point> Interpolation::findNaturalNeighbor(Point oP,Scene* sc){
    /*Centroid instertion*/
    Scene modifiedScene = *sc;
    Scene& originScene = (*sc);
    std::vector<Point> points;
    points.push_back(oP);
    std::vector<Point> neighbors;

    /*Insert new vertex*/
    int nb=modifiedScene.getVertices().size();
    Vertex_handle vertex = modifiedScene.insert_vertex(oP, 0, nb);
    if (vertex != Vertex_handle()) modifiedScene.getVertices().push_back(vertex);

    /*update centroids
     * compute modified centroids
     * update_positions();
    */
    //TO DO: add modified centroids
    unsigned j = 0;
    for (unsigned i = 0; i < nb+1; ++i)
    {
        //update modifier point
        Vertex_handle vm = modifiedScene.getVertices()[i];
        if (vm->is_hidden()) continue;
        Point pm = vm->compute_centroid();
        points.push_back(pm);
        pm = modifiedScene.getDomain().clamp(pm); //go to boundaries
        vm->set_position(pm);

        //compare original voronoi
        Vertex_handle vo = originScene.getVertices()[i];
        if (i==nb) continue;
        Point po = vo->compute_centroid();
        if (po != pm)
            neighbors.push_back(pm);
    }

    modifiedScene.update_triangulation();

    /*compute weights*/
    return neighbors;
}
