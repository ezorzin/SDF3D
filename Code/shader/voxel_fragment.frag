/// @file

#version 460 core

uniform mat4 V_mat;                                                             // View matrix.
uniform mat4 P_mat;                                                             // Projection matrix.
uniform float AR;                                                               // Framebuffer aspect ratio.

in  vec4 color;                                                                 // Voxel color.
in  vec2 quad;                                                                  // Voxel quad.

out vec4 fragment_color;                                                        // Fragment color.
 
// CONSTANTS:
#define PI 3.1415925359f                                                        // PI.
#define TWO_PI 6.2831852f                                                       // 2*PI.
#define MAX_STEPS 100                                                           // Maximum ray marching steps.
#define MAX_DISTANCE 100.0f                                                     // Maximum marching distance.
#define EPSILON 0.01f                                                           // Minimum surface distance (ray marching epsilon).
#define INF 1.0f/0.0f                                                           // Infinity.
 
//////////////////////////////////////////////////////////////////////////////////
/////////////// OBJECT'S SIGNED DISTANCE FIELDS IMPLICIT FUNCTIONS ///////////////
//////////////////////////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// SCENE RENDERING ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void main()
{
  int   i;                                                                      // Step index.
  float d;                                                                      // Marching distance.
  vec3  ray;                                                                    // Marching ray.
  float sdf = INF;                                                                    // Signed distance field.
  float sdf_L = INF;
  float sdf_R = INF;
  vec3  ro;                                                                     // Ray origin (vantage point).
  vec3  rd;                                                                     // Ray direction.
  vec3  lp;                                                                     // Light position.
  vec3  ld;                                                                     // Light direction.
  float li;                                                                     // Light intensity.
  float ambient;                                                                // Ambient light intenstity.
  float diffuse;                                                                // Diffuse light intensity.
  float specular;                                                               // Specular light intensity.
  vec3  dx = vec3(EPSILON, 0.0f, 0.0f);                                         // x-direction increment.
  vec3  dy = vec3(0.0f, EPSILON, 0.0f);                                         // y-direction increment.
  vec3  dz = vec3(0.0f, 0.0f, EPSILON);                                         // z-direction increment.
  float nx = INF;                                                               // Scene SDF's gradient x-component.
  float ny = INF;                                                               // Scene SDF's gradient y-component.
  float nz = INF;                                                               // Scene SDF's gradient z-component.
  vec3  n;                                                                      // Scene's normal vector.

  // INITIALIZING SCENE:
  lp = vec3(5.0f, 5.0f, 0.0f);                                                  // Setting light position...
  ro = vec3(0.0f, 1.0f, 0.0f);                                                  // Setting ray origin (vantage point)...
  rd = normalize(vec3(quad.x*AR, quad.y, 1.0f));                                // Computing ray direction...
  d = 0.0f;                                                                     // Resetting marching distance...
  i = 0;                                                                        // Resetting step index... 

  // RAY MARCHING:
  while ((d < MAX_DISTANCE) && (sdf > EPSILON) && (i < MAX_STEPS))
  {
    ray = ro + d*rd;                                                            // Computing marching ray...      
    sdf = INF;                                                                  // Setting SDF to infinity...
    sdf = min(sdf, SDF_plane(ray));                                             // Computing signed distance field...
    sdf = min(sdf, SDF_sphere(ray));                                            // Computing signed distance field...
    d  += sdf;                                                                  // Updating marching distance...
    i++;
  }

  ray = ro + d*rd;                                                              // Computing final marching ray...

  // COMPUTING NORMAL VECTOR:
  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(ray - dx));                                      // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(ray - dx));                                     // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(ray + dx));                                      // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(ray + dx));                                     // Computing signed distance field (right limit)...
  
  nx = sdf_R - sdf_L;                                                           // Computing gradient (x-component)...      
  
  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(ray - dy));                                      // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(ray - dy));                                     // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(ray + dy));                                      // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(ray + dy));                                     // Computing signed distance field (right limit)...
  
  ny = sdf_R - sdf_L;                                                           // Computing gradient (y-component)... 

  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(ray - dz));                                      // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(ray - dz));                                     // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(ray + dz));                                      // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(ray + dz));                                     // Computing signed distance field (right limit)...
  
  nz = sdf_R - sdf_L;                                                           // Computing gradient (z-component)... 
  
  n = normalize(vec3(nx, ny, nz));                                              // Computing normal vector...
  
  // COMPUTING LIGHTNING:
  ambient = 0.1f;
  ld = normalize(lp - ray);                                                     // Computing light direction...
  diffuse = clamp(dot(n, ld), 0.0f, 1.0f);                                      // Computing diffuse light intensity...
  li = ambient + diffuse;

  fragment_color = vec4(li, li, li, 1.0f);                                      // Setting output color...
}