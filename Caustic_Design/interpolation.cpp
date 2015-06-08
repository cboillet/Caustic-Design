#include "interpolation.h"
#include "random.h"
#include "voronoi_creation.h"

Interpolation::Interpolation(Scene* sc, Scene* tsc, Scene* csc, MainWindow* w):m_scene(sc),target_scene(tsc),compute_scene(csc), win(w){
    if (!sc->getDomain().is_valid()) return;
    double dx = sc->getDomain().get_dx();
    double dy = sc->getDomain().get_dy();
    int i;
    for (i=0; i<400; i++)
    {
        double x = random_double(-dx, dx);
        double y = random_double(-dy, dy);
        Xo.push_back( Point(x, y) );
    }
}

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
void Interpolation::runInterpolation(){
   std::vector<Vertex_handle> neighbors = findNaturalNeighbor(Xo[0]);
   std::vector<std::pair<Vertex_handle, FT> > vertices_weight = computeWeights(neighbors,Xo[0]);
  // std::vector<Point> computeLightOnDistribution = computeXr(vertices_weight);
}

std::vector<Vertex_handle> Interpolation::findNaturalNeighbor(Point oP){
    unsigned int i;
    bool newnei;
    int size;
    std::vector<Point> points;
    std::vector<Vertex_handle> neighbors;
    Vertex_handle neighborVertex;
    std::vector<Point> polygonc;
    std::vector<Point> polygono;

    /*Insert new vertex/oP as centroids as in the computational Scene*/
    std::vector<Vertex_handle>& cs_vertex = compute_scene->getVertices();
    std::cout << "compute scene vertices= " << cs_vertex.data() << std::endl;
    size = compute_scene->getVertices().size();
    std::cout << "vertices size before insertion =" << size << std::endl;
    for (i = 0; i < size; ++i)
    {
       Vertex_handle vi = m_scene->getVertices()[i];
       Vertex_handle vtest = compute_scene->getVertices()[i];

       if (vi->is_hidden()) continue;
       Point ci = vi->compute_centroid();
       points.push_back(ci);
    }
    points.push_back(oP);
    Vertex_handle vertex = compute_scene->insert_vertex(points[i], 0.0, size);
    cs_vertex.push_back(vertex);
    compute_scene->update_positions(points);
    std::vector<FT> weights(points.size(), 0.0);
    compute_scene->construct_triangulation(points, weights);
    neighbors = compute_scene->find_neighbors(compute_scene->getVertices()[size]);
    return neighbors;
}

std::vector<std::pair<Vertex_handle, FT> > Interpolation::computeWeights(std::vector<Vertex_handle> neighbors, Point oP){
    unsigned int i;
    unsigned int mscIndex;
    unsigned int cscIndex;
    std::vector<std::pair<Vertex_handle, FT> > vertices_weight;
    std::pair<Vertex_handle, FT> p;
    Vertex_handle vn;
    FT areaOnFace;
    FT weight;
    FT temp;

    if (neighbors.size() == 0) {
        std::cout << "weight vector empty" << std::endl;
        return vertices_weight;
    }

    /*Compute the ratio of area occupied by the new cell on the older cell*/
    for (i=0; i<neighbors.size(); ++i){
        /*compute area of the overlapsing cell*/
        vn = neighbors[i];
        mscIndex = m_scene->findIndexVertice(vn);
        cscIndex = compute_scene->findIndexVertice(vn);
        areaOnFace = m_scene->getVertices()[mscIndex]->compute_area() - compute_scene->getVertices()[cscIndex]->compute_area();
        /*compute ratio*/
        temp = m_scene->getVertices()[mscIndex]->compute_area();
        weight = areaOnFace/temp;
        //p.first = vn;
        //p.second = weight;
        //vertices_weight.push_back(std::make_pair(vn,weight));
    }

    return vertices_weight;
}

std::vector<Point> Interpolation::computeXr(std::vector<std::pair<Vertex_handle, FT> > vertices_weight){


}
