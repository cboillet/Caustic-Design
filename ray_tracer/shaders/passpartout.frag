#version 330

in vec3 pos;
uniform vec2 renderSize;
uniform vec2 screenSize;
uniform float opacity;

out vec4 color;

void main(void)
{
	float ratio = renderSize.y*screenSize.x/(renderSize.x*screenSize.y)*2.0;
    float relativeHeight = ratio/2.0;
    if(abs(pos.y)> relativeHeight){
		color = vec4(0.0, 0.0, 0.0,opacity);
    }else{
		color = vec4(0.0, 0.0, 0.0, 0.0);
    }
}
