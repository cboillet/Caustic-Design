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

    if(cmdline_args.count() == 4){
        filename = cmdline_args.at(1);
        output_file = cmdline_args.at(2);
        npoints = cmdline_args.at(3).toUInt();
    }else{
        std::cerr << "Usage:";
        std::cerr << "\tvoronoi-creation in_file out_file nsites" << std::endl;
        std::cerr << "\t\tin_file \t- image file that is used as source irradiance" << std::endl;
        std::cerr << "\t\tout_file \t- .dat-file to write centroid points to" << std::endl;
        std::cerr << "\t\tnsites \t\t- amount of points/sites created from image" << std::endl << std::endl;
        return 1;
    }

    VoronoiCreator vc = VoronoiCreator();
    vc.load_image(filename);
    vc.init_points(npoints);
    vc.apply_lloyd_optimization();
    vc.write_dat_file(output_file);

    return 0;

    //a.exit();

    //return a.exec();
}
