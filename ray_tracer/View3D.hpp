#ifndef VIEW3D_H
#define VIEW3D_H

#include "Model.hpp"

#include <glm/glm.hpp>

namespace View3D{
    void display(void);
    void reshape(int width, int height);

    void loadShader();

    void resetCamera();
    void updateCamera();

    extern glm::vec2 screen;

    /**
    * Glut Callback functions
    */
    void keyPressed(unsigned char key, int x, int y);

    void mouseDragged(int x, int y);

    void mousePressed(int button, int state, int x, int y);
}

#endif // VIEW3D_H
