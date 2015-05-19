#include "Ray.hpp"

#include <glm/glm.hpp>


using namespace glm;

Ray::Ray()
{
}

Ray::Ray(const Ray & ray){
    this->origin = vec3(ray.origin);
    this->dir = vec3(ray.dir);

}

Ray::Ray(vec3 origin, vec3 dir){
    this->origin = origin;
    this->dir = dir;
}


void Ray::getOrigin(vec3 &origin){
    origin = this->origin;

}

void Ray::getDirection(vec3 &dir){
    dir = this->dir;
}
