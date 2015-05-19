#version 330

in vec4 position;
in vec3 normal;

out vec3 norm;
out vec4 pos;

uniform mat4 modelViewProjectionMatrix;
uniform mat4 modelViewMatrix;
uniform mat3 normalMatrix;

void main(void)
{
	gl_Position = modelViewProjectionMatrix * position;

        norm = normalMatrix * normal;
        pos = modelViewMatrix * position;
}
