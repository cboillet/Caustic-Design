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
      zPosition = 45;
      xCenter = 0;
      drawAxis = true;
      drawNormals = true;
      drawDesiredNormals = true;
      drawDesiredRayDirections = false;
      drawDesiredRays = false;
      drawTestRay = false;
      highlightIndex = 14;
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

void Renderer::wheelEvent(QWheelEvent *event)
{
    zPosition -= float(event->delta())/200.f;
    updateCamera();
}


void Renderer::updateCamera()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (GLfloat)this->width()/(GLfloat)this->height(), 0.1f, 1000.0f);
    gluLookAt(xCenter, 0, zPosition, // eye (where camera is at)
              xCenter, 0, 0, // center (where to look at)
              0, 1, 0  // up-vector
              );
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_DEPTH_TEST);
    glLoadIdentity();

    update();
}

void Renderer::sceneUpdate()
{
    if(model.getLightRayPositions().empty())
        xCenter = 0.0f;
    else
        xCenter = 0.0 * model.getFocalLength();

    updateCamera();
}

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
    updateCamera();
}

void ModelRendering::paintGL(){
    //glCullFace(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    if(model.meshes.empty()) return;

    Mesh meshToDraw = model.meshes[0];


    // apply rotation around middle between mesh and receiver
    glTranslatef(xCenter, 0, 0);
    glRotatef(x_rotate+current_x_rotate, 1.0, 0.0, 0.0);
    glRotatef(y_rotate+current_y_rotate, 0.0, 1.0, 0.0);
    glTranslatef(-xCenter, 0, 0);

    if(drawAxis)
        paintAxis();

    paintMesh(meshToDraw);

    if(drawNormals)
        paintNormals(meshToDraw);

    if(drawDesiredNormals)
        paintDesiredNormals();

    if(drawDesiredRays)
        paintDesiredRays();

    if(drawDesiredRayDirections)
        paintDesiredRayDirections();

    paintVertexHighlight();

    paintNeighbors();

    if(drawTestRay)
        paintTestRay();

    paintReceiver();


    glFlush();
}


void ModelRendering::paintAxis()
{
    // draw axis
    glBegin(GL_LINES);
        // x
        glColor3f(0,0,1);
        glVertex3f(0.f, 0.f, 0.f);
        glVertex3f(2.0f * model.surfaceSize, 0.f, 0.f);

        // y
        glColor3f(0, 1, 0);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 2.0f * model.surfaceSize, 0);

        // z
        glColor3f(1, 0, 0);
        glVertex3f(0, 0, 0);
        glVertex3f(0, 0, 2.0f * model.surfaceSize);
    glEnd();
}


void ModelRendering::paintMesh(Mesh mesh)
{
    Vertex v;
    glBegin(GL_TRIANGLES);

    for (uint i=0; i<mesh.indices.size(); i++)
    {
        v = mesh.vertices[mesh.indices[i].x];
        glColor3f(0.0f,0.0f,1.0f);
        glVertex3f(v.Position.x, v.Position.y, v.Position.z);

        v = mesh.vertices[mesh.indices[i].y];
        glColor3f(0.0f,1.0f,0.0f);
        glVertex3f(v.Position.x, v.Position.y, v.Position.z);

        v = mesh.vertices[mesh.indices[i].z];
        glColor3f(1.0f,0.0f,0.0f);
        glVertex3f(v.Position.x, v.Position.y, v.Position.z);
    }
    glEnd();

}


void ModelRendering::paintReceiver()
{
    glPointSize(3.0f);
    glBegin(GL_POINTS);
    glColor3f(0, 1, 1);
    std::vector<glm::vec3> light_pos = model.getLightRayPositions();
    for (uint i=0; i<light_pos.size(); i++)
    {
        glm::vec3 p = light_pos[i];
        glColor3f(0,0,0);
        glVertex3f(p.x, p.y, p.z);
    }
    glEnd();

    if(!light_pos.empty())
    {
        float surfaceSize = model.surfaceSize;
        float f = model.getFocalLength() + model.meshes[0].getMaxX();
        glBegin(GL_QUADS);
            glColor3f(1,1,1);
            glVertex3f(f, -surfaceSize, -surfaceSize);
            glVertex3f(f, -surfaceSize, surfaceSize);
            glVertex3f(f, surfaceSize, surfaceSize);
            glVertex3f(f, surfaceSize, -surfaceSize);
        glEnd();
    }
}

void ModelRendering::paintNormals(Mesh mesh)
{
    glBegin(GL_LINES);
    glColor3f(1,1,1);
    for (uint i=0; i<mesh.vertices.size(); i++)
    {
        glm::vec3 pos = mesh.vertices[i].Position;
        glm::vec3 end = pos + mesh.vertices[i].Normal;
        glVertex3f(pos.x, pos.y, pos.z);
        glVertex3f(end.x, end.y, end.z);
    }
    glEnd();
}

void ModelRendering::paintDesiredNormals()
{
    glBegin(GL_LINES);
    glColor3f(0,1,1);
    for (uint i=0; i<model.desiredNormals.size(); i++)
    {
        glm::vec3 pos = model.meshes[0].faceVertices[i]->Position;
        glm::vec3 end = model.desiredNormals[i] + pos;
        glVertex3f(pos.x, pos.y, pos.z);
        glVertex3f(end.x, end.y, end.z);
    }

    glEnd();
}

void ModelRendering::paintDesiredRayDirections()
{
    glBegin(GL_LINES);

    glColor3f(0,0,1);
    for (uint i=0; i<model.screenDirections.size(); i++)
    {
        glm::vec3 pos = model.meshes[0].faceVertices[i]->Position;
        glm::vec3 end = model.screenDirections[i] + pos;
        glVertex3f(pos.x, pos.y, pos.z);
        glVertex3f(end.x, end.y, end.z);
    }

    glEnd();
}

void ModelRendering::paintDesiredRays()
{
    if(model.meshes[0].faceVertices.empty() || model.receiverLightPositions.empty())
        return;

    glBegin(GL_LINES);
    glColor3f(0, 1, 0);

    Mesh m = model.meshes[0];
    for (uint i=0; i<model.meshes[0].faceVertices.size(); i++)
    {
        glVertex3f(m.faceVertices[i]->Position.x, m.faceVertices[i]->Position.y, m.faceVertices[i]->Position.z);
        glVertex3f(model.receiverLightPositions[i].x, model.receiverLightPositions[i].y, model.receiverLightPositions[i].z);
    }
    glEnd();
}

void ModelRendering::paintVertexHighlight()
{
    if(true) return;

    if(highlightIndex != -1 && !model.meshes[0].vertices.empty())
    {

        Mesh mesh = model.meshes[0];
        Vertex* v = &mesh.vertices[highlightIndex];

        //std::vector<int> neighbors = mesh.getNeighborsIndex(v, false);

        glPointSize(10.0f);
        glBegin(GL_POINTS);
        glColor3f(1,1,0);

        glVertex3f(v->Position.x, v->Position.y, v->Position.z);

        glColor3f(0,1,1);
        for(uint i=0; i<neighbors.size(); i++)
        {
            int index = neighbors[i];
            glm::vec3 pos = mesh.vertices[index].Position;
            glVertex3f(pos.x, pos.y, pos.z);
        }

        glEnd();
    }
}

void ModelRendering::paintNeighbors()
{
    int vertexIndex = 14;

    for(uint i=0; i<neighborMapping.size(); i+=2)
    {
        float color[3] = {0, 0, 0};
        int index = (i%6)/2;
        color[index] = 1;

        glColor3f(color[0], color[1], color[2]);
        Vertex* v1 = &model.meshes[0].vertices[vertexIndex];
        Vertex* v2 = &model.meshes[0].vertices[neighbors[neighborMapping[i]]];
        Vertex* v3 = &model.meshes[0].vertices[neighbors[neighborMapping[i+1]]];

        glm::vec3 xoffset = glm::vec3(0.1, 0, 0) * float(i);

        glBegin(GL_LINES);
        // first edge
        glm::vec3 pos = v1->Position + xoffset;
        glVertex3f(pos.x, pos.y, pos.z);
        pos = v2->Position + xoffset;
        glVertex3f(pos.x, pos.y, pos.z);

        // second edge
        glVertex3f(pos.x, pos.y, pos.z);
        pos = v3->Position + xoffset;
        glVertex3f(pos.x, pos.y, pos.z);

        // third edge
        glVertex3f(pos.x, pos.y, pos.z);
        pos = v1->Position + xoffset;
        glVertex3f(pos.x, pos.y, pos.z);
        glEnd();

        // highlight starting-point (= end-point)
        glBegin(GL_POINTS);
        glVertex3f(pos.x, pos.y, pos.z);

        // and 2nd point
        pos = v2->Position + xoffset;
        glVertex3f(pos.x, pos.y, pos.z);
        glEnd();

    }

}

void ModelRendering::paintTestRay()
{
    glm::vec3 startPos = glm::vec3(-10, 0, 0);
    glBegin(GL_LINES);
    glColor3f(1, 0.33, 0.9);

    glVertex3f(startPos.x, startPos.y, startPos.z);
    glVertex3f(testRayRedirect.x, testRayRedirect.y, testRayRedirect.z);

    glVertex3f(testRayRedirect.x, testRayRedirect.y, testRayRedirect.z);
    glVertex3f(testRayEndPoint.x, testRayEndPoint.y, testRayEndPoint.z);

    glEnd();

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
