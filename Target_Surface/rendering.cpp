#include "global.h"

#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>

// Qt
#include <QtGui>
#include <QDialog>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QClipboard>

#include "rendering.h"
#include "SurfaceModel.h"
#include <iostream>

Renderer::Renderer(QWidget *parent):QGLWidget(parent){
    y_rotate = 0.0f;
    mouse_is_down = false;
}

Renderer::Renderer(int framesPerSecond, QWidget *parent , char *name):QGLWidget(parent){
    setWindowTitle(QString::fromUtf8(name));
      if(framesPerSecond == 0)
          t_Timer = NULL;
      else
      {
          int seconde = 1000; // 1 seconde = 1000 ms
          int timerInterval = seconde / framesPerSecond;
          t_Timer = new QTimer(this);
          connect(t_Timer, SIGNAL(timeout()), this, SLOT(timeOutSlot()));
          t_Timer->start( timerInterval );
      }

      y_rotate = 0.0f;
      mouse_is_down = false;
}

void Renderer::printVersion(){
    std::cout<<"gpu"<<glGetString(GL_RENDERER)<<"opengl version:"<< glGetString(GL_VERSION) <<std::endl;
}

void Renderer::mousePressEvent(QMouseEvent * evt)
{
    mouse_down = new QMouseEvent(*evt);
    mouse_is_down = true;
}

void Renderer::keyPressEvent(QKeyEvent *keyEvent)
{
    switch(keyEvent->key())
    {
        case Qt::Key_Escape:
            close();
            break;
    case Qt::Key_R:
            y_rotate = 0.0f;
            break;
    }
}

void Renderer::mouseReleaseEvent(QMouseEvent *)
{
    mouse_is_down = false;
    delete mouse_down;
    mouse_down = NULL;
    y_rotate += current_y_rotate;
    current_y_rotate = 0;
    x_rotate += current_x_rotate;
    current_x_rotate = 0;
}

void Renderer::mouseMoveEvent(QMouseEvent * evt)
{
    if(mouse_is_down)
    {
        current_y_rotate = 0.5f* float(evt->pos().x() - mouse_down->pos().x());
        current_x_rotate = 0.5f* float(evt->pos().y() - mouse_down->pos().y());
        update();
    }
}

void Renderer::timeOutSlot(){}



/*Model Rendering*/

ModelRendering::ModelRendering(QWidget *parent):Renderer(60,parent,"Model Rendering")
{
    //setModel();
    //model.meshes[0].generateTriangles();

}

void ModelRendering::initializeGL(){
    initializeGLFunctions();
    glShadeModel(GL_SMOOTH);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    //glEnable(GL_CULL_FACE);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

}


void ModelRendering::resizeGL(int width, int height){
    if(height == 0)
        height = 1;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (GLfloat)width/(GLfloat)height, 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_DEPTH_TEST);
    glLoadIdentity();
}
void ModelRendering::paintGL(){
    //glCullFace(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    if(model.meshes.empty()) return;

    Mesh meshToDraw = model.meshes[0];
    //drawMesh(meshToDraw);
    Vertex v;
    int i=0;

    glTranslatef(0.0f, 0.0f, -6.0f);
    glRotatef(x_rotate+current_x_rotate, 1.0, 0.0, 0.0);
    glRotatef(y_rotate+current_y_rotate, 0.0, 1.0, 0.0);


    // draw axis
    glBegin(GL_LINES);
        // x
        glColor3f(0,0,1);
        glVertex3f(0.f, 0.f, 0.f);
        glVertex3f(2.0f, 0.f, 0.f);

        // y
        glColor3f(0, 1, 0);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 2, 0);

        // z
        glColor3f(1, 0, 0);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 0, 2);
    glEnd();

    glBegin(GL_TRIANGLES);
    while (i<meshToDraw.vertices.size()){                   //la ligne coupe le triangle et est partiellement cachée
        v = meshToDraw.vertices[i];
        glColor3f(0.0f,0.0f,1.0f);
        glVertex3f(v.Position.x, v.Position.y, v.Position.z);

        v = meshToDraw.vertices[i+1];
        glColor3f(0.0f,1.0f,0.0f);
        glVertex3f(v.Position.x, v.Position.y, v.Position.z);

        v = meshToDraw.vertices[i+2];
        glColor3f(1.0f,0.0f,0.0f);
        glVertex3f(v.Position.x, v.Position.y, v.Position.z);

        i+=3;
     }

    glEnd();
    glFlush();

//    Mesh meshToDraw = model.meshes[0];
//    glGenVertexArrays(1, &meshToDraw.VAO);
//    glGenBuffers(1, &meshToDraw.VBO);
//    glGenBuffers(1, &meshToDraw.EBO);

//    glBindVertexArray(meshToDraw.VAO);
//    glBindBuffer(GL_ARRAY_BUFFER, meshToDraw.VBO);

//    glBufferData(GL_ARRAY_BUFFER, meshToDraw.vertices.size() * sizeof(Vertex),
//                 &meshToDraw.vertices[0], GL_STATIC_DRAW);

//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, meshToDraw.EBO);
//    glBufferData(GL_ELEMENT_ARRAY_BUFFER, meshToDraw.indices.size() * sizeof(GLuint),
//                 &meshToDraw.indices[0], GL_STATIC_DRAW);

//    // Vertex Positions
//    glEnableVertexAttribArray(0);
//    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
//                         (GLvoid*)0);
//    // Vertex Normals
//    glEnableVertexAttribArray(1);
//    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
//                         (GLvoid*)offsetof(Vertex, Normal));
//    // Vertex Texture Coords
//    glEnableVertexAttribArray(2);
//    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex),
//                         (GLvoid*)offsetof(Vertex, TexCoords));

//    glBindVertexArray(0);

//    glBindVertexArray(meshToDraw.VAO);
//    glDrawElements(GL_TRIANGLES, meshToDraw.indices.size(), GL_UNSIGNED_INT, 0);
//    glBindVertexArray(0);


//    glBegin(GL_TRIANGLES);
//        glVertex3f(0.0f, 1.0f, 0.0f);
//        glVertex3f(-1.0f, -1.0f, 0.0f);
//        glVertex3f(1.0f, -1.0f, 0.0f);
//    glEnd();

//    glTranslatef(3.0f, 0.0f, -6.0f);


}


void ModelRendering::setModel(){
    QString modelToLoad = QFileDialog::getOpenFileName(this, tr("Open 3D model to load"), ".");
    if(modelToLoad.isEmpty()) return;
    model.loadModel(modelToLoad.toStdString());
    //setUpMesh(model.meshes[0]);
    update();
    model.printAllVertices();
    initializeGLFunctions();
    //glewInit(); //not any more necessary with qglfunction

}
void ModelRendering::setUpMesh(Mesh meshToDraw){
    makeCurrent();
    std::cout << "setupMesh function" << std::endl;
    glGenVertexArrays(1, &meshToDraw.VAO);
    glGenBuffers(1, &meshToDraw.VBO);
    glGenBuffers(1, &meshToDraw.EBO);

    glBindVertexArray(meshToDraw.VAO);
    glBindBuffer(GL_ARRAY_BUFFER, meshToDraw.VBO);

    glBufferData(GL_ARRAY_BUFFER, meshToDraw.vertices.size() * sizeof(Vertex),
                 &meshToDraw.vertices[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, meshToDraw.EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, meshToDraw.indices.size() * sizeof(GLuint),
                 &meshToDraw.indices[0], GL_STATIC_DRAW);

    // Vertex Positions
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                         (GLvoid*)0);
    // Vertex Normals
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                         (GLvoid*)offsetof(Vertex, Normal));
    // Vertex Texture Coords
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                         (GLvoid*)offsetof(Vertex, TexCoords));

    glBindVertexArray(0);

}

void ModelRendering::drawMesh(Mesh meshToDraw){
    makeCurrent();
    glBindVertexArray(meshToDraw.VAO);
    glDrawElements(GL_TRIANGLES, meshToDraw.indices.size(), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

}
