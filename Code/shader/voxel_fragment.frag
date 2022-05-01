/// @file

#version 460 core

uniform mat4 V_mat;                                                                                 // View matrix.
uniform mat4 P_mat;                                                                                 // Projection matrix.
uniform float AR;                                                                                   // Framebuffer aspect ratio.

in  vec4 color;                                                                                     // Voxel color.
in  vec2 quad;                                                                                      // Voxel quad.

out vec4 fragment_color;                                                                            // Fragment color.
 
// CONSTANTS:
#define PI 3.1415925359f                                                                            // PI.
#define TWO_PI 6.2831852f                                                                           // 2*PI.
#define MAX_STEPS 100                                                                               // Maximum ray marching steps.
#define MAX_DISTANCE 100.0f                                                                         // Maximum marching distance.
#define EPSILON 0.01f                                                                               // Minimum surface distance (ray marching epsilon).
#define INF 1.0f/0.0f                                                                               // Infinity.
#define DX vec3(EPSILON, 0.0f, 0.0f);                                                               // x-direction increment.
#define DY vec3(0.0f, EPSILON, 0.0f);                                                               // y-direction increment.
#define DZ vec3(0.0f, 0.0f, EPSILON);                                                               // z-direction increment.

// CAMERA:
struct Camera
{
  vec3  pos;                                                                                        // Camera position.
  float fov;                                                                                        // Camera field of view (degrees).
};

// LIGHT:
struct Light
{
  vec3  pos;                                                                                        // Light position.
  vec3  col;                                                                                        // Light color.
  float amb;                                                                                        // Light ambient intensity.
  float dif;                                                                                        // Light diffused intensity.          
  float ref;                                                                                        // Light reflected intensity.
  vec3  amb_col;                                                                                    // Light ambient color.
  vec3  dif_col;                                                                                    // Light diffused color.
  vec3  ref_col;                                                                                    // Light reflected color.
};

// MATERIAL:
struct Material
{
  vec3  amb_col;                                                                                    // Material ambient color.
  vec3  dif_col;                                                                                    // Material diffused color.
  vec3  ref_col;                                                                                    // Material reflected color.
  float sh;                                                                                         // Material shininess.
};

//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// OBJECT'S SIGNED DISTANCE FIELDS IMPLICIT FUNCTIONS /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
float sphereSDF(vec3 position)
{
  float x = 0.0f;                                                                                   // Sphere x-coordinate.
  float y = 0.3f;                                                                                   // Sphere y-coordinate.
  float z = 0.0f;                                                                                   // Sphere z-coordinate.
  float r = 0.2f;                                                                                   // Sphere radius.
  vec3  s = vec3(x, y, z);                                                                          // Sphere position.
  float sdf = length(position - s) - r;                                                             // Signed distance field.
 
  return sdf;                                                                                       // Returning signed distance field...
}

float planeSDF(vec3 position)
{
  float sdf = position.y;                                                                           // Signed distance field.
 
  return sdf;                                                                                       // Returning signed distance field...
}

float sceneSDF(vec3 position)
{
  float sdf = INF;

  sdf = min(sdf, planeSDF(position));                                                               // Computing signed distance field...
  sdf = min(sdf, sphereSDF(position));                                                              // Computing signed distance field...

  return sdf;                                                                                       // Returning signed distance field...
}

float raymarch(vec3 position, vec3 direction)
{
  float distance = 0.0f;
  float sdf = INF;
  float i = 0;
  vec3  ray;

  while (
          (distance < MAX_DISTANCE) &&                                                              // Checking marching distance...
          (sdf > EPSILON) &&                                                                        // Checking signed distance field...
          (i < MAX_STEPS)                                                                           // Checking step index...
        )
  {
    ray = position + distance*direction;                                                            // Computing marching ray...      
    distance += sdf;                                                                                // Updating marching distance...
    i++;                                                                                            // Updating step index...
  }

  return distance;                                                                                  // Returning final distance...
}

vec3 normal(vec3 position)
{
  float nx;                                                                                         // Scene SDF's gradient x-component.
  float ny;                                                                                         // Scene SDF's gradient y-component.
  float nz;                                                                                         // Scene SDF's gradient z-component.
  float sdf_L;                                                                                      // Scene SDF (left limit).
  float sdf_R;                                                                                      // Scene SDF (right limit).

  sdf_L = sceneSDF(position - DX);                                                                  // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DX);                                                                  // Computing signed distance field (left limit)...
  nx = sdf_R - sdf_L;                                                                               // Computing gradient (x-component, centered derivative)...

  sdf_L = sceneSDF(position - DY);                                                                  // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DY);                                                                  // Computing signed distance field (left limit)...
  ny = sdf_R - sdf_L;                                                                               // Computing gradient (x-component, centered derivative)...

  sdf_L = sceneSDF(position - DZ);                                                                  // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DZ);                                                                  // Computing signed distance field (left limit)...
  nz = sdf_R - sdf_L;                                                                               // Computing gradient (x-component, centered derivative)...   
  
  return (normalize(vec3(nx, ny, nz)));                                                               // Returnig normal vector...
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// SCENE RENDERING ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
void main()
{
  Camera    camera;
  Light     light;
  vec3      ray;
  vec3      P;
  vec3      N;
  float     d;
  Material  M;

  // INITIALIZING RAY MARCHING:
  camera.fov = 60;                                                                                  // Setting camera field of view...
  camera.pos = vec3(0.0f, 0.1f, 0.6f);                                                              // Setting camera position...
  camera.pos = (inverse(V_mat)*vec4(camera.pos, 1.0f)).xyz;                                         // Applying arcball to camera position...

  light.pos = vec3(5.0f, 5.0f, 0.0f);                                                               // Setting light position...
  light.col = vec3(0.7f, 0.7f, 0.7f);                                                               // Setting light color...
  light.amb = 0.1f;                                                                                 // Setting light ambient intensity...

  ray = normalize(vec3(quad.x*AR, quad.y, -2.0f/tan(camera.fov*PI/360.0f)));                        // Computing position on canvas (ray intersection on quad)...
  ray = normalize(inverse(V_mat)*vec4(ray, 0.0f)).xyz;                                              // Applying arcball to camera direction...
  
  M.amb_col   = vec3(0.0f, 0.2f, 0.8f);                                                             // Setting material ambient color...
  M.dif_col   = vec3(0.0f, 0.2f, 0.8f);                                                             // Setting material diffuse color...
  M.ref_col   = vec3(0.5f, 0.5f, 0.5f);                                                             // Setting material specular color...
  M.sh        = 32.0f;                                                                              // Setting material shininess...

  // COMPUTING RAY MARCHING:
  d = raymarch(camera.pos, ray);                                                                    // Computing scene distance...
  P = camera.pos + d*ray;                                                                           // Computing scene position...
  N = normal(P);                                                                                    // Computing scene normal...
  
  // COMPUTING LIGHTNING:
  ray = normalize(light.pos - ray);                                                                 // Computing light direction...
  light.dif = clamp(dot(N, ray), 0.0f, 1.0f);                                                       // Computing light diffused intensity...
  light.ref = pow(max(dot(normalize(camera.pos - ray), reflect(-ray, N)), 0.0f), M.sh);             // Computing light reflected intensity...
  light.amb_col = light.amb*light.col;                                                              // Computing light ambient color...
  light.dif_col = light.col*light.dif;                                                              // Computing light diffused color...
  light.ref_col = light.col*light.ref*M.ref_col;                                                    // Computing reflected light color...

  /*
  // COMPUTING RAY MARCHING FOR SHADOWS:
  while (
          (distance < MAX_DISTANCE) &&                                                              // Checking matching distance...
          (sdf > EPSILON) &&                                                                        // Checking signed distance field...
          (i < MAX_STEPS)                                                                           // Checking step index...
        )
  {
    ray_march = (p + n*EPSILON*2.0f) + distance*normalize(light_position - p);                      // Computing marching ray...      
    sdf = INF;                                                                                      // Setting SDF to infinity...
    sdf = min(sdf, SDF_plane(ray_march));                                                           // Computing signed distance field...
    sdf = min(sdf, SDF_sphere(ray_march));                                                          // Computing signed distance field...
    distance += sdf;                                                                                // Updating marching distance...
    i++;                                                                                            // Updating step index...
  }

  ray_march = (p + n*EPSILON*2.0f) + distance*normalize(light_position - p);                       // Computing final marching ray...

  if(distance < length(ray_march))
  {
    light_diffused_intensity *= 0.1f;
  }
  */

  //light_diffuse_color = light_color*(light_diffused_intensity*material_diffused_color);             // Computing diffuse light color...
  

  fragment_color = vec4(light.amb_col + light.dif_col + light.ref_col, 1.0f);                       // Setting output color...
}