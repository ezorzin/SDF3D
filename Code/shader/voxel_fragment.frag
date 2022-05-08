/// @file

#version 460 core

uniform mat4 V_mat;                                                                                 // View matrix.
uniform mat4 P_mat;                                                                                 // Projection matrix.
uniform float AR;                                                                                   // Framebuffer aspect ratio.

in  vec2 quad;                                                                                      // Voxel quad.

layout(std430, binding = 0) buffer voxel_color
{
  vec4 color_SSBO[1024*768];                                                                         // Voxel color SSBO.
};

out vec4 fragment_color;                                                                            // Fragment color.

void main()
{
  fragment_color = color_SSBO[int((quad.x + 1.0f)*(1024 - 1)/2.0f) + 1024*int((quad.y + 1.0f)*(768 - 1)/2.0f)];       // Setting fragment color...
}