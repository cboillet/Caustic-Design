#include "interpolation.h"
#include "random.h"
#include "voronoi_creation.h"
#include "config.h"
#include <math.h>


std::mutex mtx;

Interpolation::Interpolation(Scene* sc, Scene* tsc, Scene* csc, int sitesAmount, MainWindow* w):m_scene(sc),source_scene(tsc),compute_scene(csc), win(w){
    if (!sc->getDomain().is_valid()) return;
    if (!tsc->getDomain().is_valid()) return;

    std::cout << "nbpoints" << source_scene->getLightPointsSource().size() << std::endl;
    Xo=source_scene->getLightPointsSource();
}


void Interpolation::runInterpolation(QString imageFile, QString datFile)
{
    Timer::start_timer(compute_scene->m_timer, COLOR_BLUE, "Interpolation(mutli-thread)");

    unsigned int n = std::thread::hardware_concurrency();
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

    source_scene->getLightPointsTarget().resize(source_scene->getLightPointsSource().size());
    m_scene->getLightPointsTarget().resize(source_scene->getLightPointsSource().size());

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


void multiThreadInterpolation(Interpolation* inter, uint id, uint from, uint to)
{

    uint range = to - from;

    for (uint i = from; i<to; i++)
    {
        std::cout << "thread-state: " << (i-from) << "/" << range << std::endl;
        findNeighborr(inter->Xo[i], inter, id, i);
    }
}

void findNeighborr(Point oP, Interpolation* inter, uint id, int index)
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

    int* points = new int[coords.size()];

    FT x = 0, y = 0;

    for(uint i=0; i<coords.size(); i++)
    {

        std::map<Point, Vertex_handle>::iterator it = inter->map.find(Point(coords[i].first.x(), coords[i].first.y()));
        Vertex_handle vh = it->second;
        int index = vh->get_index();
        Point centroid = inter->m_scene->m_vertices[index]->compute_centroid();
        x += centroid.x() * coords[i].second;
        y += centroid.y() * coords[i].second;

    }

    x /= norm;
    y /= norm;

    Point p = Point(x,y);

    mtx.lock();
    inter->m_scene->getLightPointsTarget()[index] = p;
    inter->source_scene->getLightPointsTarget()[index] = p;
    mtx.unlock();

    delete[] points;
}

