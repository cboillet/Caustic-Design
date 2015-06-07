#include "interpolation.h"
#include "random.h"
#include "voronoi_creation.h"

<<<<<<< HEAD
Interpolation::Interpolation(Scene* sc, Scene* tsc, Scene* csc):m_scene(sc),target_scene(tsc),compute_scene(csc){
=======
Interpolation::Interpolation(Scene* sc){
>>>>>>> 13c14360646e9cad8f42118e91dc7f290592b60b
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

<<<<<<< HEAD
bool Interpolation::prepareData(){
  /*  std::vector<Vertex_handle> compute_vertices;
    std::vector<FT> compute_weights;
    std::vector<Point> compute_points;

    //ensure scene are available
    if(!m_scene) return false;
    if(!compute_scene) return false;
    std::cout << "scenes available.. ";

    // --- retrieve points, weights, vertices
    source_points.clear();
    std::vector<FT> scene_weights = std::vector<FT>();
    m_scene->collect_sites(source_points, scene_weights);

    compute_points.clear();
    compute_weights.clear();
    compute_scene->collect_sites(compute_points, compute_weights);

    source_vertices = m_scene->getVertices();
    compute_vertices = compute_scene->getVertices();

    // --- ensure they are of same dimension
    //if(target_points.size() != source_points.size()) return false;
    std::cout << "same point amount.. ";
    //if(target_weights.size() != scene_weights.size()) return false;
    std::cout << "same weight amount.. ";
    if(compute_vertices.size() != source_vertices.size())
    {
        std::cout << "error.. target_vertices.size = " << target_vertices.size() << " != " << source_vertices.size() << " = source.vertices.size";
        return false;
    }
    std::cout << "same vertex amount.. ";

    std::cout << std::endl;
    // --- no issue found
    return true;
    */
}

std::vector<Vertex_handle> Interpolation::compareCell(Vertex_handle vc){
   /* Edge_circulator ecirc = compute_scene->getRT().incident_edges(vi);
    Edge_circulator eend  = ecirc;
    Vertex_handle neighborVertex;
    std::vector<Vertex_handle> neighbors;
    GAL_For_all(ecirc, eend)
    {
        Edge edge = *ecirc;
        neighborVertex = get_opposite(edge);
        std::cout << "we have a neighbor here" << std::endl;
        neighbors.push_back(neighborVertex);
    }


    std::vector<Point> polygonc;
    std::vector<Point> polygono;
    if (!vc->is_hidden()){
        bool ok = compute_scene->getRT().pre_build_polygon(vc, vc->dual().points());
        compute_scene->getRT().build_polygon(vc, polygonc);
    }

    for (int i = 0; i <= m_scene->getVertices().size(); ++i){
        compute_scene->getRT().build_polygon(vo, polygono);
        if(polygono == polygonc) {
            std::cout << "not a neigbor" << std::endl;
            return false;
        }
        polygono.clear();
    }

    std::cout << "we have a neighbor here" << std::endl;
    return true;
    */

}

=======
>>>>>>> 13c14360646e9cad8f42118e91dc7f290592b60b
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
}

<<<<<<< HEAD
std::vector<Vertex_handle> Interpolation::findNaturalNeighbor(Point oP){
    unsigned int i;
    std::vector<Point> points;
    std::vector<Vertex_handle> neighbors;
    Vertex_handle neighborVertex;
    bool newnei;
    /*Insert new vertex/oP as centroids as in the computational Scene*/
    //test
    if(!compute_scene){
        std::cerr << "target scene not available!" << std::endl;
    }
    std::vector<Vertex_handle> cs_vertex = compute_scene->getVertices();
    std::cout << "compute scene vertices= " << cs_vertex.data() << std::endl;

    Vertex_handle testc=cs_vertex[3];
    for (i = 0; i < m_scene->getVertices().size(); ++i)
    {
       Vertex_handle vi = m_scene->getVertices()[i];
       Vertex_handle vtest = compute_scene->getVertices()[i];

       if (vi->is_hidden()) continue;
       Point ci = vi->compute_centroid();
       points.push_back(ci);
    }
    points.push_back(oP);
    Vertex_handle vertex = compute_scene->insert_vertex(points[i], 0.0, cs_vertex.size());
    cs_vertex.push_back(vertex);
    compute_scene->update_positions(points);
    std::vector<FT> weights(points.size(), 0.0);
    compute_scene->construct_triangulation(points, weights);

    //need not centroids but edge !

    /*reconstruct the cell polygon and check the modified polygon*/
    /*
      for (i = 0; i <= m_scene->getVertices().size(); ++i){
       if(points[i] != oP )

           newnei = this->compareCell(m_scene->getVertices()[i], compute_scene->getVertices()[i]);
       }
    */
    const RT& rt = compute_scene->getRT();
    Edge_circulator ecirc = rt.incident_edges(vertex); //to debug !
   /* Edge_circulator eend  = ecirc;
    CGAL_For_all(ecirc, eend)
    {
        Edge edge = *ecirc;
        neighborVertex = compute_scene->getRT().get_opposite(edge);
        std::cout << "we have a neighbor here" << std::endl;
        neighbors.push_back(neighborVertex);

    }
    */


    std::vector<Point> polygonc;
    std::vector<Point> polygono;
    if (!vertex->is_hidden()){
        bool ok = compute_scene->getRT().pre_build_polygon(vertex, vertex->dual().points());
        compute_scene->getRT().build_polygon(vertex, polygonc);
    }

    for (int i = 0; i <= m_scene->getVertices().size(); ++i){
        compute_scene->getRT().build_polygon(m_scene->getVertices()[i], polygono);
        if(polygono == polygonc) {
            std::cout << "not a neigbor" << std::endl;
        }
        polygono.clear();
    }

    std::cout << "we have a neighbor here" << std::endl;
    return neighbors;
=======

std::vector<Point> Interpolation::findNaturalNeighbor(Point oP,Scene* sc){
    unsigned int i;
    std::vector<Point> points;
//    for (i = 0; i < sc->getVertices().size(); ++i)
//    {
//        Vertex_handle vi = sc->getVertices()[i];
//        if (vi->is_hidden()) continue;
//        Point ci = vi->compute_centroid();
//        points.push_back(ci);
//    }
    //points.push_back(oP);
   // Vertex_handle vertex = modifiedScene.insert_vertex(points[i+1], 0.0, modifiedScene.getVertices().size());
   // modifiedScene.getVertices().push_back(vertex);
    //std::vector<FT> weights(points.size(), 0.0);
    //modifiedScene.construct_triangulation(points, weights);

    return points;
>>>>>>> 13c14360646e9cad8f42118e91dc7f290592b60b
}


//std::vector<Point> Interpolation::findNaturalNeighbor(Point oP,Scene* sc){
//    /*Centroid instertion*/
//    Scene modifiedScene = *sc;
//    Scene& originScene = (*sc);
//    std::vector<Point> points;
//    points.push_back(oP);
//    std::vector<Point> neighbors;

//    /*Insert new vertex*/
//    int nb=modifiedScene.getVertices().size();
//   // Vertex_handle vertex = modifiedScene.insert_vertex(oP, 0, nb);
//   // if (vertex != Vertex_handle()) modifiedScene.getVertices().push_back(vertex);
//    std::cout << "recomputing Voronoi for the copy";
//    std::cout << std::flush;
//    modifiedScene.optimize_positions_via_lloyd(true);
//    std::cout << "done" << std::endl;

//    /*update centroids
//     * compute modified centroids
//     * update_positions();
//    */
//    //TO DO: add modified centroids


//    /*
//    {
//        //update modifier point
//        Vertex_handle vm = modifiedScene.getVertices()[i];
//        if (vm->is_hidden()) continue;
//        Point pm = vm->compute_centroid();
//        points.push_back(pm);
//        pm = modifiedScene.getDomain().clamp(pm); //go to boundaries
//        vm->set_position(pm);

//        //compare original voronoi
//        Vertex_handle vo = originScene.getVertices()[i];
//        if (i==nb) continue;
//        Point po = vo->compute_centroid();
//        if (po != pm)
//            neighbors.push_back(pm);
//    }

//    modifiedScene.update_positions(points);
//    modifiedScene.update_triangulation(); //new centroid inserted
//    */
//    /*compute weights*/
//    return neighbors;
//}
