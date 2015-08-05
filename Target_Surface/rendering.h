#ifndef RENDERING_H
#define RENDERING_H
#include "global.h"
#include "SurfaceModel.h"
#include "SurfaceMesh.h"
#include <iostream>

class Renderer : public QGLWidget, protected QGLFunctions{
Q_OBJECT

public:
    Renderer(QWidget *parent);
    Renderer(int framesPerSecond, QWidget *parent, char *name);
    virtual void initializeGL(){printVersion();}
    virtual void resizeGL(int width, int height) {}
    virtual void paintGL() {}
    virtual void mouseMoveEvent(QMouseEvent * evt);
    virtual void mousePressEvent(QMouseEvent * evt);
    virtual void mouseReleaseEvent(QMouseEvent * evt);
    void printVersion();
    void keyPressEvent(QKeyEvent *keyEvent);

    float y_rotate;
    float x_rotate;
    float current_y_rotate;
    float current_x_rotate;

    Model model;

public slots:
    virtual void timeOutSlot();

private:
    QTimer *t_Timer;
    //Shader shader;
    QMouseEvent* mouse_down;
    bool mouse_is_down;

};

class ModelRendering : public Renderer{
Q_OBJECT

public:
    ModelRendering(QWidget *parent = 0);
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();

    void setModel();
    void setUpMesh(Mesh meshToDraw);
    void drawMesh(Mesh meshToDraw);
};

#endif // RENDERING_H
