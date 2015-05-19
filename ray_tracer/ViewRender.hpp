#ifndef VIEWRENDER_H
#define VIEWRENDER_H

#include "glm/glm.hpp"
#include <vector>

namespace ViewRender{
    void display(void);
    void reshape(int width, int height);
    void mouseDragged(int x, int y);
    void mousePressed(int button, int state, int x, int y);
    void keyPressed(unsigned char key, int x, int y);

    void buffer(void);
    void init();
    void updateImage();
    void addPoint(glm::vec3 pointPosition, glm::vec3 pointColor);
    void normalize();
    void removeOldPoints();
    void updateCamera();
    void resetCamera();

    extern glm::vec2 screen;
}


#endif // VIEWRENDER_H
