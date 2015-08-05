#include "mainwindow.h"
#include <QApplication>


using namespace glm;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    glewExperimental = GL_TRUE;
    glewInit();

    //ModelRendering mr;

    MainWindow w;
    //mr.show();
    w.show();
    return a.exec();
}


