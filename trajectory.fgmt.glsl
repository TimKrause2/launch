#version 320 es
precision highp float;

layout(location =0) out vec4 outColor;

uniform vec4 color;

void main(){

    // Output color = color of the texture at the specified UV
    outColor = color;
}
