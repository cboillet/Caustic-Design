#ifndef VORONOICREATOR_H
#define VORONOICREATOR_H

#include <scene.h>
#include <QtCore>
#include <QString>

class VoronoiCreator
{
public:
    VoronoiCreator();

    ~VoronoiCreator(){
        delete m_scene;
    }

    void init_points(int npoints);
    void load_image(QString filename);
    void apply_lloyd_optimization();
    void write_dat_file(QString output_file);

    QString output_file;
    Scene* m_scene;
};

#endif // VORONOICREATOR_H
