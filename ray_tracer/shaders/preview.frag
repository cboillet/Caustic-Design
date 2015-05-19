#version 330

in vec4 pos;
in vec3 norm;

out vec4 color;

vec4 blinnPhongReflection(vec3 position, vec3 normal);

void main(void)
{
        color = blinnPhongReflection(pos.xyz, normalize(norm));
}
