/// @file

#version 460 core

uniform mat4 V_mat;                                                             // View matrix.
uniform mat4 P_mat;                                                             // Projection matrix.

in  vec4 color;                                                                 // Voxel color.
in  vec2 quad;                                                                  // Voxel quad.

out vec4 fragment_color;                                                        // Fragment color.

void main(void)
{
  fragment_color = vec4(1.0f, 1.0f, 1.0f, 1.0f);                                // Setting fragment color...  
}