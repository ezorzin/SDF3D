/// @file

#version 460 core

uniform mat4 V_mat;                                                             // View matrix.
uniform mat4 P_mat;                                                             // Projection matrix.
uniform float AR;                                                               // Framebuffer aspect ratio.

in  vec4 color;                                                                 // Voxel color.
in  vec2 quad;                                                                  // Voxel quad.

out vec4 fragment_color;                                                        // Fragment color.
 
// Constants
#define PI 3.1415925359f                                                        // PI.
#define TWO_PI 6.2831852f                                                       // 2*PI.
#define MAX_STEPS 100                                                           // Maximum ray marching steps.
#define MAX_DIST 100.0f                                                         // Maximum distance.
#define SURFACE_DIST 0.01f                                                      // Surface distance.
 
float SDF_sphere(vec3 p)
{
  float x = 0.0f;                                                               // Sphere x-coordinate.
  float y = 1.0f;                                                               // Sphere y-coordinate.
  float z = 6.0f;                                                               // Sphere z-coordinate.
  float r = 1.0f;                                                               // Sphere radius.
  vec3  s = vec3(x, y, z);                                                      // Sphere position.
  float SDF = length(p - s) - r;                                                // Signed distance field.
 
  return SDF;                                                                   // Returning signed distance field...
}

float SDF_plane(vec3 p)
{
  float SDF = p.y;                                                              // Signed distance field.
 
  return SDF;                                                                   // Returning signed distance field...
}
 
float RayMarching(vec3 ro, vec3 rd) 
{
  float dO = 0.0f;                                                              // Distance origin
    
  for(int i=0;i<MAX_STEPS;i++)
  {
    vec3 p = ro + rd * dO;
    float ds = GetDist(p);                                                      // ds is Distance Scene
    dO += ds;
    
    if(dO > MAX_DIST || ds < SURFACE_DIST) break;
  }
  
  return dO;
}
 
vec3 GetNormal(vec3 p)
{ 
  float d_plane = SDF_plane(p);
  float d_sphere = SDF_sphere(p);
  float d = min(d_plane, d_sphere);

  vec2 e = vec2(0.01f, 0.0f); // Epsilon
  vec3 n = d - vec3(GetDist(p-e.xyy), GetDist(p-e.yxy), GetDist(p-e.yyx));
   
  return normalize(n);
}

float GetLight(vec3 p)
{ 
  // Light (directional diffuse)
  vec3 lightPos = vec3(5.0f, 5.0f, 0.0f); // Light Position
  vec3 l = normalize(lightPos-p); // Light Vector
  vec3 n = GetNormal(p); // Normal Vector
   
  float dif = dot(n,l); // Diffuse light
  dif = clamp(dif,0.,1.); // Clamp so it doesnt go below 0
 
  return dif;
}

void main()
{
  vec3 ro = vec3(0,1,0);                                                        // Ray Origin/ Camera
  vec3 rd = normalize(vec3(quad.x*AR,quad.y,1));                                // Ray Direction
  float d = RayMarching(ro,rd);                                                 // Distance field
  vec3 p = ro + rd * d;
  float dif = GetLight(p); // Diffuse lighting
  d*= .2;
  vec3 color = vec3(dif);
  
  // Set the output color
  fragment_color = vec4(color, 1.0);
}