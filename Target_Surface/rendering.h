#ifndef RENDERING_H
#define RENDERING_H
#include "global.h"
#include "SurfaceModel.h"
#include "SurfaceMesh.h"
#include <iostream>

class Renderer : public QGLWidget, protected QGLFunctions{
Q_OBJECT

public:
    Renderer(int framesPerSecond, QWidget *parent, char *name);
    virtual void initializeGL(){printVersion();}
    virtual void resizeGL(int width, int height) {}
    virtual void paintGL() {}
    virtual void keyPressEvent(QKeyEvent *keyEvent);
    void printVersion();


public slots:
    virtual void timeOutSlot();

private:
    QTimer *t_Timer;
    //Shader shader;

};

class ModelRendering : public Renderer{
Q_OBJECT

public:
    ModelRendering(QWidget *parent = 0);
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();

    Model model;
    void setModel();
    void setUpMesh(Mesh meshToDraw);
    void drawMesh(Mesh meshToDraw);
};

#endif // RENDERING_H
