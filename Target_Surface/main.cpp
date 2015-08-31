#include "mainwindow.h"
#include <QApplication>

using namespace glm;

int main(int argc, char *argv[])
{
    google::InitGoogleLogging("logtostderr=1");

    QApplication a(argc, argv);
    glewExperimental = GL_TRUE;
    glewInit();

    std::cout << "running with EDIR " << EDIR_WEIGHT << ", EREG " << EREG_WEIGHT << std::endl;

    //ModelRendering mr;

    MainWindow w;
    //mr.show();
    w.show();
    return a.exec();
}


