#include "interpolation.h"
#include "random.h"
#include "voronoi_creation.h"
#include "config.h"
#include <math.h>


std::mutex mtx;

Interpolation::Interpolation(Scene* sc, Scene* tsc, Scene* csc, int sitesAmount, MainWindow* w):m_scene(sc),source_scene(tsc),compute_scene(csc), win(w){
    if (!sc->getDomain().is_valid()) return;
    if (!tsc->getDomain().is_valid()) return;

//    int nbWith = floor(MESH_AMOUNT) / 3;
//    int nbHeight = MESH_AMOUNT / nbWith;
//    FT stepx = 2.0 * tsc->getDomain().get_dx() / nbWith;
//    FT stepy = 2.0 * tsc->getDomain().get_dy() / nbHeight;
//    std::cout << "dx" << tsc->getDomain().get_dy() << std::endl;
//    std::cout << "dy " << tsc->getDomain().get_dy() << std::endl;
//    std::cout << "stepx" << stepx << std::endl;
//    std::cout << "stepy" << stepy << std::endl;
    int nbpoints = 0;
    int scene_sites = sitesAmount;

//    source_scene->getLightPointsSource().clear();
//    for (unsigned i = 0; i < nbWith; ++i)
//    {
//        FT x = (i + 0.5)*stepx - tsc->getDomain().get_dx();
//        x += EPS;
//        for (unsigned j = 0; j < nbHeight; ++j)
//        {
//            FT y = (j + 0.5)*stepy - tsc->getDomain().get_dy();
//            y += EPS;
//            Xo.push_back(Point(x, y));
//            source_scene->getLightPointsSource().push_back(Point(x,y));
//            std::cout << "inserted point Xrs x" << x << std::endl;
//            std::cout << "inserted point Xrs y" << y << std::endl;
//            nbpoints++;
//        }
//    }
    std::cout << "nbpoints" << source_scene->getLightPointsSource().size() << std::endl;
    Xo=source_scene->getLightPointsSource();
    //prepare compute scene
    VoronoiCreator voronoi_creator;
}


void Interpolation::runInterpolation(QString imageFile, QString datFile)
{
    Timer::start_timer(compute_scene->m_timer, COLOR_BLUE, "Interpolation(mutli-thread)");

    unsigned int n = 3;//std::thread::hardware_concurrency();
    if(n==0) n=4;

    computeScenes = new Scene*[n];
    computeScenes[0] = compute_scene;

    uint i;
    for(i=1; i<n; i++)
    {
        computeScenes[i] = new Scene();

        computeScenes[i]->load_image(imageFile);
        computeScenes[i]->load_points(datFile);
    }

    std::cout<<"light point source size"<<source_scene->getLightPointsSource().size()<<std::endl;

    for(i=0; i < compute_scene->m_vertices.size(); i++)
    {
         Point p = compute_scene->m_vertices[i]->get_position();
         Vertex_handle v = compute_scene->m_vertices[i];
         map.insert(std::make_pair(p,v));
    }

    std::thread* ts = new std::thread[n];

    int step = source_scene->getLightPointsSource().size() / n;
    int max = source_scene->getLightPointsSource().size();
    for(uint j=0; j<n; j++)
    {
         uint from, to;
         from = j*step;
         to = (j+1)*step;

         if(to > max)
             to = max;

         if(j == (n-1)){
             if(to < max)
                 to = max;
         }

         ts[j] = std::thread(multiThreadInterpolation, this, j, from, to);
    }

    for (uint j=0; j<n; j++)
    {
        ts[j].join();
    }

    delete[] ts;

    delete[] computeScenes;
    Timer::stop_timer(compute_scene->m_timer, COLOR_BLUE);
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
    Timer::start_timer(compute_scene->m_timer, COLOR_BLUE, "Interpolation");
    std::cout << std::endl;

   int i;
   std::cout<<"light point source size"<<source_scene->getLightPointsSource().size()<<std::endl;

   for(i=0; i < compute_scene->m_vertices.size(); i++)
   {
        Point p = compute_scene->m_vertices[i]->get_position();
        Vertex_handle v = compute_scene->m_vertices[i];
        map.insert(std::make_pair(p,v));
   }

   for(i=0; i<source_scene->getLightPointsSource().size(); i++)
   {
       findNeighbor(Xo[i]);
   }

   Timer::stop_timer(compute_scene->m_timer, COLOR_BLUE);

   win->viewer->toggle_view_Xrs();
   win->viewer->toggle_view_Xr();
   win->update();

   if(true)
       return;

   int promille = -1;
   int currentPromille = -1;

   for(i=0; i<source_scene->getLightPointsSource().size(); ++i)
   {

       /*currentPromille = int(float(i)*1000.0f / float(source_scene->getLightPointsSource().size()));
       if(currentPromille > promille)
       {
           promille = currentPromille;
           std::cout << "\rcurrent promille: " << promille;
           std::cout.flush();
       }*/

       std::cout << "\rlight point " << i << "/" << source_scene->getLightPointsSource().size();
       std::cout.flush();

       findNeighbor(Xo[i]);
   }
   std::cout << std::endl;

   Timer::stop_timer(compute_scene->m_timer, COLOR_BLUE);


   // std::vector<Point> computeLightOnDistribution = computeXr(vertices_weight);
   /*next: 1. insert weight
           2. compute weight for all Xo -> compute surface delete point/reset as before insertion
           3. Xo normal distribution
           4. UI for output
    */
   win->viewer->toggle_view_Xrs();
   win->viewer->toggle_view_Xr();
   win->update();
}

void multiThreadInterpolation(Interpolation* inter, uint id, uint from, uint to)
{

    uint range = to - from;

    for (uint i = from; i<to; i++)
    {
        std::cout << "thread-state: " << (i-from) << "/" << range << std::endl;
        findNeighborr(inter->Xo[i], inter, id);
    }
}

void findNeighborr(Point oP, Interpolation* inter, uint id)
{
    RT rt = inter->computeScenes[id]->m_rt;

    Weighted_point wp(oP, 0);
    Point_coordinate_vector coords;
    CGAL::Triple<
            std::back_insert_iterator<Point_coordinate_vector>,
            FT, bool> result =
            CGAL::regular_neighbor_coordinates_2(
                rt,
                wp,
                std::back_inserter(coords));

    if(!result.third)
    {
        std::cout << "The coordinate computation was not successful."
                  << std::endl;
        std::cout << "The point (" <<wp.point() << ") lies outside the convex hull."
                  << std::endl;
    }

    FT  norm = result.second;
      //std::cout << "Coordinate computation successful." << std::endl;
      //std::cout << "Normalization factor: " <<norm << std::endl;
      //std::cout << "done" << std::endl;

    int* points = new int[coords.size()];

    FT x = 0, y = 0;
    //std::cout << "in-point " << wp << std::endl;
    for(uint i=0; i<coords.size(); i++)
    {
        //std::cout << "coords: " << coords[i].first << ", weight: " << coords[i].second << std::endl;

        std::map<Point, Vertex_handle>::iterator it = inter->map.find(Point(coords[i].first.x(), coords[i].first.y()));
        Vertex_handle vh = it->second;
        int index = vh->get_index();
        //int index = compute_scene->findIndexByPoint(coords[i].first);
        Point centroid = inter->m_scene->m_vertices[index]->compute_centroid();
        x += centroid.x() * coords[i].second;
        y += centroid.y() * coords[i].second;

        //points[i] = compute_scene->findIndexByPoint(coords[i].first);
        //std::cout << "point index: " << points[i] << std::endl;
    }
    //std::cout << std::endl;

    x /= norm;
    y /= norm;

    Point p = Point(x,y);

    mtx.lock();
    inter->m_scene->getLightPointsTarget().push_back(p);
    inter->source_scene->getLightPointsTarget().push_back(p);
    mtx.unlock();

    delete[] points;
}

void Interpolation::findNeighbor(Point oP)
{
    RT rt = compute_scene->m_rt;

    Weighted_point wp(oP, 0);
    Point_coordinate_vector coords;
    CGAL::Triple<
            std::back_insert_iterator<Point_coordinate_vector>,
            FT, bool> result =
            CGAL::regular_neighbor_coordinates_2(
                rt,
                wp,
                std::back_inserter(coords));

    if(!result.third)
    {
        std::cout << "The coordinate computation was not successful."
                  << std::endl;
        std::cout << "The point (" <<wp.point() << ") lies outside the convex hull."
                  << std::endl;
    }

    FT  norm = result.second;
      //std::cout << "Coordinate computation successful." << std::endl;
      //std::cout << "Normalization factor: " <<norm << std::endl;
      //std::cout << "done" << std::endl;

    int* points = new int[coords.size()];

    FT x = 0, y = 0;
    //std::cout << "in-point " << wp << std::endl;
    for(uint i=0; i<coords.size(); i++)
    {
        //std::cout << "coords: " << coords[i].first << ", weight: " << coords[i].second << std::endl;

        std::map<Point, Vertex_handle>::iterator it = map.find(Point(coords[i].first.x(), coords[i].first.y()));
        Vertex_handle vh = it->second;
        int index = vh->get_index();
        //int index = compute_scene->findIndexByPoint(coords[i].first);
        Point centroid = m_scene->m_vertices[index]->compute_centroid();
        x += centroid.x() * coords[i].second;
        y += centroid.y() * coords[i].second;

        //points[i] = compute_scene->findIndexByPoint(coords[i].first);
        //std::cout << "point index: " << points[i] << std::endl;
    }
    //std::cout << std::endl;

    x /= norm;
    y /= norm;

    Point p = Point(x,y);

    m_scene->getLightPointsTarget().push_back(p);
    source_scene->getLightPointsTarget().push_back(p);

    delete[] points;
}


void Interpolation::findNaturalNeighbor(Point oP){
    unsigned int i;
    bool newnei;
    int size;
    std::vector<Point> points;
    std::vector<Vertex_handle> neighbors;
    Vertex_handle neighborVertex;

    /*Insert new vertex/oP as centroids as in the computational Scene*/
    std::vector<Vertex_handle>& cs_vertex = compute_scene->getVertices();
    size = cs_vertex.size();
    //size = compute_scene->m_vertices.size();


    for (i = 0; i < size; ++i)
    {

       Vertex_handle vi = source_scene->getVertices()[i];
       if (vi->is_hidden()) continue;
       Point ci = vi->compute_centroid();
       points.push_back(ci);

    }

    points.push_back(oP);

    Vertex_handle vertex = compute_scene->insert_vertex(points[i], 0.0, size);
    compute_scene->m_vertices.push_back(vertex);
    cs_vertex.push_back(vertex);

    //compute_scene->update_triangulation_values();

    compute_scene->update_positions(points);
    std::vector<FT> weights(points.size(), 0.0);
    compute_scene->construct_triangulation(points, weights);

    neighbors = compute_scene->find_neighbors(compute_scene->getVertices()[size]);

    computeWeights(neighbors, oP);

    /*reset the computational scene*/
    //compute_scene->delete_vertex(vertex);
    points.clear();
    neighbors.clear();
    cs_vertex.pop_back();
    //compute_scene->m_vertices.pop_back();
    return;
}


void Interpolation::computeWeights(std::vector<Vertex_handle> neighbors, Point oP){
    unsigned int i;
    unsigned int mscIndex, cscIndex;
    std::vector<std::pair<Vertex_handle, FT> > vertices_weight;
    std::pair<Vertex_handle, FT> p;
    Vertex_handle vn;
    double x = 0;
    double y = 0;
//    FT areaOnFace;
//    FT weight;
//    FT temp;

    if (neighbors.size() == 0) {
        std::cout << "weight vector empty" << std::endl;
        return;
    }

    /*Compute the ratio of area occupied by the new cell on the older cell*/
    for (i=0; i<neighbors.size(); ++i){
        /*compute area of the overlapsing cell*/
        vn = neighbors[i];
        mscIndex = source_scene->findIndexVerticeBySite(vn);
        cscIndex = compute_scene->findIndexVerticeByCentroid(vn);
        std::cout<<"index source scene"<<mscIndex<<std::endl;
        std::cout<<"index compute scene"<<cscIndex<<std::endl;

        // FT area_source_scene = source_scene->getVertices()[mscIndex]->compute_area();
       // FT area_c_scene = compute_scene->getVertices()[cscIndex]->compute_area();

        FT area_source_scene = source_scene->getVertices()[cscIndex]->compute_area();
        FT area_c_scene = compute_scene->getVertices()[cscIndex]->compute_area();
        std::cout<<"area_source_scene: "<<area_source_scene<<std::endl;
        std::cout<<"area_c_scene: "<<area_c_scene<<std::endl;

        FT areaOnFace= area_source_scene-area_c_scene;
        /*compute ratio*/
        FT temp = source_scene->getVertices()[mscIndex]->compute_area();
        FT weight = areaOnFace/temp;
        std::cout<<"weight"<<weight<<std::endl;
        /*retrieve on m_scene*/
        x = x + (m_scene->getVertices()[mscIndex]->compute_centroid().x())*weight;
        y = y + (m_scene->getVertices()[mscIndex]->compute_centroid().y())*weight;
    }
    std::cout<<"x for Xr"<<x<<std::endl;
    std::cout<<"y for Xr"<<y<<std::endl;
    m_scene->getLightPointsTarget().push_back(Point(x,y));
    m_scene->getLightPointsTarget().push_back(Point(x,y));
    source_scene->getLightPointsTarget().push_back(Point(x,y));
    source_scene->getLightPointsTarget().push_back(Point(x,y));

    return;
}

