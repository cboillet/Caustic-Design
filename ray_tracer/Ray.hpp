#ifndef RAY_HPP
#define RAY_HPP

#include <glm/glm.hpp>

class Ray
{
public:
    Ray();
    Ray(const Ray & ray);
    Ray(glm::vec3 origin, glm::vec3 dir);
    void getOrigin(glm::vec3 &orign);
    void getDirection(glm::vec3 &dir);

protected:
    glm::vec3 origin;
    glm::vec3 dir;
};

#endif // RAY_HPP
