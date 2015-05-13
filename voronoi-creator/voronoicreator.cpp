#include "voronoicreator.h"
#include <sstream>
#include <iostream>
#include <ostream>

VoronoiCreator::VoronoiCreator():
    m_scene(new Scene)
{
}

void VoronoiCreator::init_points(int npoints){
    std::cout << "initializing " << npoints << " points...";
    std::cout << std::flush;
    m_scene->generate_random_sites_based_on_image(npoints);
    std::cout << "done" << std::endl;
}

void VoronoiCreator::load_image(QString filename){
    std::cout << "opening " << filename.toUtf8().constData() << "... ";
    std::cout << std::flush;
    m_scene->load_image(filename);
    std::cout << "done" << std::endl;
}

void VoronoiCreator::apply_lloyd_optimization(){
    std::cout << "running lloyd optimization... ";
    std::cout << std::flush;
    m_scene->optimize_positions_via_lloyd(true);
    std::cout << "done" << std::endl;
}

void VoronoiCreator::write_dat_file(QString output_file){
    std::cout << "writing file to " << output_file.toUtf8().constData() << "... ";
    std::cout << std::flush;
    m_scene->save_points(output_file);
    std::cout << "done" << std::endl;
}
