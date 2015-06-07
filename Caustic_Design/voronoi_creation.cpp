#include "voronoi_creation.h"
#include <sstream>
#include <iostream>
#include <ostream>

void VoronoiCreator::init_points(int npoints,Scene* sc){
    std::cout << "initializing " << npoints << " points...";
    std::cout << std::flush;
    sc->generate_random_sites_based_on_image(npoints);
    std::cout << "done" << std::endl;
}

void VoronoiCreator::apply_lloyd_optimization(Scene* sc){
    std::cout << "running lloyd optimization... ";
    std::cout << std::flush;
    sc->optimize_positions_via_lloyd(true);
    std::cout << "done" << std::endl;
}


