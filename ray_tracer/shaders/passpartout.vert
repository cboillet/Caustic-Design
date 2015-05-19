#version 330

in vec3 position;
out vec3 pos;

void main(void)
{
    gl_Position = vec4(position, 1.0);
    pos = position;
}

