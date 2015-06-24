#include "voronoi_creation.h"
#include "config.h"
#include <sstream>
#include <iostream>
#include <ostream>

bool VoronoiCreator::generate_voronoi(Scene *sc, unsigned npoints, double epsilon)
{
    // --- initialize the voronoi diagram
    init_points(npoints, sc);

    FT threshold = sc->compute_position_threshold(epsilon);
    bool success = false;
    unsigned iter = 0;
    unsigned max_iter = 200;

    std::cout << "optimizing voronoi diagram via lloyd ..";
    // optimize until the movement is below a certain threshold
    while(iter++ < LLYOD_STEPS){
        std::cout << iter << " ..";

        // norm is the norm of the position gradient
        // ( a metric for how far the centroid is away from the site after the optimzation step )
        FT norm = sc->optimize_positions_via_lloyd(true);

        if(norm < threshold)
        {
            success = true;
            break;
        }
    }

    std::cout << std::endl;

    if(success)
        std::cout << "voronoi diagram created with " << iter << " optimization steps" << std::endl;
    else
        std::cerr << "voronoi could not be optimized" << std::endl;

    return success;
}

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


