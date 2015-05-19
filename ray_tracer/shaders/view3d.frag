#version 330

// Interpolated values from the vertex shaders
//in vec2 UV;
in vec4 out_color;

// Ouput data
out vec4 color;

//uniform sampler2D textureSampler;

void main(){
        color = out_color;
        //color = texture2D (textureSampler, UV);
}
