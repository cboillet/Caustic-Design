#include <iostream>
#include <QCoreApplication>
#include <QStringList>
#include <QString>

#include <voronoicreator.h>

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    a.setApplicationName("Caustic Design - Irradiance to Voronoi");

    QStringList cmdline_args = QCoreApplication::arguments();

    QString filename = QString();
    QString output_file = QString();
    unsigned int npoints = 0;
    unsigned int iteration_loops = 0;

    if(cmdline_args.count() == 5){
        filename = cmdline_args.at(1);
        output_file = cmdline_args.at(2);
        npoints = cmdline_args.at(3).toUInt();
        iteration_loops = cmdline_args.at(4).toUInt();
    }else{
        std::cerr << "Usage:";
        std::cerr << "\tvoronoi-creation in_file out_file n_sites n_lloyd" << std::endl;
        std::cerr << "\t\tin_file \t- image file that is used as source irradiance" << std::endl;
        std::cerr << "\t\tout_file \t- .dat-file to write centroid points to" << std::endl;
        std::cerr << "\t\tn_sites\t\t- amount of points/sites created from image" << std::endl;
        std::cerr << "\t\tn_lloyd\t\t- amount of iterations to run lloyd optimization" << std::endl << std::endl;
        return 1;
    }

    VoronoiCreator vc = VoronoiCreator();
    vc.load_image(filename);
    vc.init_points(npoints);
    for (uint i=0; i<iteration_loops; i++){
        std::cout << "(" << (i+1) << "/" << iteration_loops << "): ";
        vc.apply_lloyd_optimization();
    }
    vc.write_dat_file(output_file);

    return 0;

    //a.exit();

    //return a.exec();
}
