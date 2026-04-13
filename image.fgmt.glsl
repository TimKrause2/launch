#version 320 es
precision highp float;

// Interpolated values from the vertex shaders
in vec2 UV;

// Values that stay constant for the whole mesh.
uniform sampler2D myTextureSampler;
uniform float alpha;

layout(location =0) out vec4 outColor;

void main(){

        // Output color = color of the texture at the specified UV
        vec4 texel = texture( myTextureSampler, UV );
        texel.a *= alpha;
        outColor = texel;
}
