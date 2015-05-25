#version 330

//in vec2 texcoord;
in vec4 position;
in vec4 color;

//out vec2 UV;
out vec4 out_color;

uniform mat4 viewProjectionMatrix;

void main(){

    // Output position of the vertex, in clip space : MVP * position
    gl_Position =  viewProjectionMatrix * position;
    out_color = color;
    // UV of the vertex. No special space for this one.
    //UV = texcoord;
}
