/// @file
// CONSTANTS:
#define PI 3.1415925359f                                                                            // PI.
#define TWO_PI 6.2831852f                                                                           // 2*PI.
#define MAX_STEPS 100                                                                               // Maximum ray marching steps.
#define MAX_DISTANCE 100.0f                                                                         // Maximum marching distance.
#define EPSILON 0.01f                                                                               // Minimum surface distance (ray marching epsilon).
#define INF 1.0f/0.0f                                                                               // Infinity.
#define DX vec3(EPSILON, 0.0f, 0.0f)                                                                // x-direction increment.
#define DY vec3(0.0f, EPSILON, 0.0f)                                                                // y-direction increment.
#define DZ vec3(0.0f, 0.0f, EPSILON)                                                                // z-direction increment.

// CAMERA:
struct Camera
{
  float3  pos;                                                                                        // Camera position.
  float fov;                                                                                        // Camera field of view (degrees).
};

// LIGHT:
struct Light
{
  float3  pos;                                                                                        // Light position.
  float3  col;                                                                                        // Light color.
  float amb;                                                                                        // Light ambient intensity.
  float dif;                                                                                        // Light diffused intensity.
  float ref;                                                                                        // Light reflected intensity.
};

// MATERIAL:
struct Material
{
  float3  amb;                                                                                        // Material ambient color.
  float3  dif;                                                                                        // Material diffused color.
  float3  ref;                                                                                        // Material reflected color.
  float shn;                                                                                        // Material shininess.
};

//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// OBJECT'S SIGNED DISTANCE FIELDS IMPLICIT FUNCTIONS /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
float sphereSDF(float3 position)
{
  float x = 0.0f;                                                                                   // Sphere x-coordinate.
  float y = 0.4f;                                                                                   // Sphere y-coordinate.
  float z = 0.0f;                                                                                   // Sphere z-coordinate.
  float r = 0.2f;                                                                                   // Sphere radius.
  float3  s = (float3)(x, y, z);                                                                          // Sphere position.
  float sdf = length(position - s) - r;                                                             // Signed distance field.
 
  return sdf;                                                                                       // Returning signed distance field...
}

float planeSDF(float3 position)
{
  float sdf = position.y;                                                                           // Signed distance field.
 
  return sdf;                                                                                       // Returning signed distance field...
}

float sceneSDF(float3 position)
{
  float sdf = INF;                                                                                  // Signed distance field.

  sdf = min(sdf, planeSDF(position));                                                               // Computing signed distance field...
  sdf = min(sdf, sphereSDF(position));                                                              // Computing signed distance field...

  return sdf;                                                                                       // Returning signed distance field...
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// RAY MARCHING ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
float raymarch(float3 position, float3 direction)
{
  float distance = 0.0f;                                                                            // Marching distance.
  float sdf;                                                                                        // Signed distance field.
  int   i;                                                                                          // Step index.
  float3  ray;                                                                                        // Marching ray.

  for (i = 0; i < MAX_STEPS; i++)
  {
    ray = position + distance*direction;                                                            // Computing marching ray...    
    sdf = sceneSDF(ray);                                                                            // Computing scene distance field...
    distance += sdf;                                                                                // Updating marching distance...

    if(distance > MAX_DISTANCE || sdf < EPSILON) break;                                             // Checking numeric precision constraints...
  }

  return distance;                                                                                  // Returning final marching distance...
}

float shadow(float3 position, float3 direction, float k)
{
  float distance = 0.0f;                                                                            // Marching distance.
  float sdf = INF;                                                                                  // Previous step signed distance field.
  float sdf_new;                                                                                    // Current step signed distance field.
  int   i;                                                                                          // Step index.
  float3  ray;                                                                                        // Marching ray.
  float shadow = 1.0f;                                                                              // Shadow intensity.
  float intersection;                                                                               // Previous to current SDF on-ray intersection point.
  float d_est;                                                                                      // Estimated closest distance.

  for (i = 0; i < MAX_STEPS; i++)
  {
    ray = position + distance*direction;                                                            // Computing marching ray...   
    sdf_new = sceneSDF(ray);                                                                        // Computing scene distance field... 
    intersection = (i == 0) ? 0.0f : sdf_new*sdf_new/(2.0f*sdf);                                    // Computing on-ray intersection point...
    d_est = sqrt(sdf_new*sdf_new - intersection*intersection);                                      // Computing estimated closest distance...
    shadow = min(shadow, k*d_est/max(0.0f, distance - intersection));                               // Computing shadow intensity...
    sdf = sdf_new;                                                                                  // Backing up signed distance field...
    distance += sdf_new;                                                                            // Updating marching distance...
    
    if(distance > MAX_DISTANCE || shadow < EPSILON) break;                                          // Checking numeric precision constraints...
  }

  shadow = clamp(shadow, 0.0f, 1.0f);                                                               // Clamping shadow intensity...

  return shadow;                                                                                    // Returning final shadow intensity...
}

float3 normal(float3 position)
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
  
  return normalize((float3)(nx, ny, nz));                                                               // Returning normal vector...
}

__kernel void thekernel(__global float4*    color,                              // Color.
                        __global float4*    position                            // Position.
                        )
{ 
  struct Camera    camera;                                                                                 // Camera.
  struct Light     light;                                                                                  // Light.
  struct Material  M;                                                                                      // Material.
  float3      ray;                                                                                    // Marching ray.
  float     d;                                                                                      // Marching distance.
  float3      view;                                                                                   // View direction.
  float3      incident;                                                                               // Incident light direction.
  float3      reflected;                                                                              // Reflected light direction.
  float3      halfway;                                                                                // Halfway light direction (Blinn-Phong).
  float3      ambient;                                                                                // Ambient light color.
  float3      diffusion;                                                                              // Diffusion light color.
  float3      reflection;                                                                             // Reflection light color.
  float3      P;                                                                                      // Scene position.
  float3      N;                                                                                      // Scene normal direction.
  
  // INITIALIZING RAY MARCHING:
  camera.fov = 60;                                                                                  // Setting camera field of view...
  camera.pos = vec3(0.0f, 0.2f, 2.0f);                                                              // Setting camera position...
  //camera.pos = (inverse(V_mat)*vec4(camera.pos, 1.0f)).xyz;                                         // Applying arcball to camera position...

  light.pos = vec3(5.0f, 5.0f, 0.0f);                                                               // Setting light position...
  light.col = vec3(0.7f, 0.7f, 0.7f);                                                               // Setting light color...
  light.amb = 0.1f;                                                                                 // Setting light ambient intensity...

  M.amb = vec3(0.0f, 0.2f, 0.8f);                                                                   // Setting material ambient color...
  M.dif = vec3(0.0f, 0.2f, 0.8f);                                                                   // Setting material diffuse color...
  M.ref = vec3(0.5f, 0.5f, 0.5f);                                                                   // Setting material specular color...
  M.shn   = 12.0f;                                                                                  // Setting material shininess...

  ray = normalize(vec3(quad.x*AR, quad.y, -2.0f/tan(camera.fov*PI/360.0f)));                        // Computing position on canvas (ray intersection on quad)...
  ray = normalize((inverse(V_mat)*vec4(ray, 0.0f)).xyz);                                            // Applying arcball to camera direction...
  
  // COMPUTING RAY MARCHING:
  d = raymarch(camera.pos, ray);                                                                    // Computing scene distance...
  P = camera.pos + d*ray;                                                                           // Computing scene position...
  N = normal(P);                                                                                    // Computing scene normal...
  
  // COMPUTING LIGHTNING:
  view = normalize(camera.pos - P);                                                                 // Computing view direction...
  incident = normalize(light.pos - P);                                                              // Computing incident light direction...
  reflected = reflect(-incident, N);                                                                // Computing reflected light direction...
  halfway = normalize(incident + view);                                                             // Coputing halfway vector (Blinn-Phong)...
  light.ref = pow(max(dot(N, halfway), 0.0f), M.shn);                                               // Computing light reflection intensity
  light.dif = clamp(dot(N, incident), 0.0f, 1.0f)*shadow(P + N*2.0f*EPSILON, incident, 10.0f);      // Computing light diffusion intensity... 
  ambient = light.amb*M.amb;                                                                        // Computing light ambient color...
  diffusion = light.dif*M.dif;                                                                      // Computing light diffused color...
  reflection = light.ref*M.ref;                                                                     // Computing light reflected color...

  fragment_color = vec4(ambient + diffusion + reflection, 1.0f);                                    // Setting output color...
}
