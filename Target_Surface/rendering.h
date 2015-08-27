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
    void wheelEvent(QWheelEvent * event);
    void printVersion();
    void keyPressEvent(QKeyEvent *keyEvent);
    void updateCamera();
    void sceneUpdate();
    void toggleDrawNormals(){drawNormals=!drawNormals;}
    void toggleDrawDesiredNormals(){drawDesiredNormals=!drawDesiredNormals;}
    void toggleDrawAxis(){drawAxis=!drawAxis;}
    void toggleDrawDesiredRays(){drawDesiredRays = !drawDesiredRays;}
    void toggleDrawDesiredRayDirections(){drawDesiredRayDirections=!drawDesiredRayDirections;}
    void setNeigbors(vector<int> neighbors, std::vector<int> neighborMapping){
        this->neighbors = neighbors;
        this->neighborMapping = neighborMapping;
    }

    float y_rotate;
    float x_rotate;
    float current_y_rotate;
    float current_x_rotate;

    float zPosition;
    float xCenter;

    Model model;

public slots:
    virtual void timeOutSlot();

protected:
    bool drawDesiredRays;
    bool drawDesiredRayDirections;
    bool drawNormals;
    bool drawDesiredNormals;
    bool drawAxis;

    std::vector<int> neighbors;
    std::vector<int> neighborMapping;

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
    void paintAxis();
    void paintMesh(Mesh mesh);
    void paintReceiver();
    void paintNormals(Mesh mesh);
    void paintDesiredNormals();
    void drawHighlights();
    void paintDesiredRays();
    void paintDesiredRayDirections();
    void paintNeighbors();
    void paintEdgeVertices();
    void paintVertex();

    void setModel();
    void setUpMesh(Mesh meshToDraw);
    void drawMesh(Mesh meshToDraw);
};

#endif // RENDERING_H
