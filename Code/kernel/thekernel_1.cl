/// @file
// CONSTANTS:
#define PI 3.1415925359f                                                                            // PI.
#define TWO_PI 6.2831852f                                                                           // 2*PI.
#define MAX_STEPS 100                                                                               // Maximum ray marching steps.
#define MAX_DISTANCE 100.0f                                                                         // Maximum marching distance.
#define EPSILON 0.01f                                                                               // Minimum surface distance (ray marching epsilon).
#define INF 1.0f/0.0f                                                                               // Infinity.
#define DX (float3)(EPSILON, 0.0f, 0.0f)                                                            // x-direction increment.
#define DY (float3)(0.0f, EPSILON, 0.0f)                                                            // y-direction increment.
#define DZ (float3)(0.0f, 0.0f, EPSILON)                                                            // z-direction increment.

// CAMERA:
struct Camera
{
  float3 pos;                                                                                       // Camera position.
  float  fov;                                                                                       // Camera field of view (degrees).
};

// LIGHT:
struct Light
{
  float3  pos;                                                                                      // Light position.
  float4  col;                                                                                      // Light color.
};

// MATERIAL:
struct Material
{
  float4  amb;                                                                                      // Material ambient color.
  float4  dif;                                                                                      // Material diffused color.
  float4  ref;                                                                                      // Material reflected color.
  float4  shn;                                                                                      // Material shininess.
};

// FRAGMENT:
struct Fragment
{
  float4  amb;
  float4  dif;
  float4  ref;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// OBJECT'S SIGNED DISTANCE FIELDS IMPLICIT FUNCTIONS /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
float sphereSDF(float3 position)
{
  float x = 0.0f;                                                                                   // Sphere x-coordinate.
  float y = 0.0f;                                                                                   // Sphere y-coordinate.
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
  
  return normalize((float3)(nx, ny, nz));                                                           // Returning normal vector...
}

float3 reflect(float3 I, float3 N)
{
  float3 reflect;

  reflect = I - 2.0f*dot(N, I)*N;

  return reflect;
}

float4 mul(float4* M, float4 V)
{
  float4 A;

  A.x = dot(M[0], V);
  A.y = dot(M[1], V);
  A.z = dot(M[2], V);
  A.w = dot(M[3], V);

  return A;
}

__kernel void thekernel(__global float4*    V_mat,                                                  // View matrix [4x4].
                        __global float4*    canvas,                                                 // Canvas [W, H, AR, FOV].
                        __global float4*    light_position,                                         // Light position [x, y, z, 1.0f].
                        __global float4*    light_color,                                            // Light color [r, g, b, a].
                        __global float4*    material,                                               // Material.
                        __global float4*    fragment_color                                          // Fragment color.
                        )
{ 
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////// INDICES //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  uint            i = get_global_id(0);                                                             // Global index i [#].
  uint            j = get_global_id(1);                                                             // Global index j [#].

  float4          V[4];                                                                             // View matrix.

  float           W = canvas[0].x;                                                                  // Window width [px].
  float           H = canvas[0].y;                                                                  // Window height [px].
  float           AR = canvas[0].z;                                                                 // Window aspect ratio [].
  float           FOV = canvas[0].w;                                                                // FOV [degrees].
  float           x = 2.0f*(i/W) - 1.0f;                                                            // Canvas x-coordinate [-1.0f...+1.0f].
  float           y = 2.0f*(j/H) - 1.0f;                                                            // Canvas y-coordinate [-1.0f...+1.0f].

  struct Camera   camera;                                                                           // Camera.
  struct Light    light;                                                                            // Light.
  struct Material M;                                                                                // Material.
  struct Fragment fragment;                                                                         // Fragment.

  float3          ray;                                                                              // Marching ray.
  float           d;                                                                                // Marching distance.
  float3          view;                                                                             // View direction.
  float3          incident;                                                                         // Incident light direction.
  float3          reflected;                                                                        // Reflected light direction.
  float3          halfway;                                                                          // Halfway light direction (Blinn-Phong).
  float           ambient;                                                                          // Ambient light intensity.
  float           diffusion;                                                                        // Diffusion light intensity.
  float           reflection;                                                                       // Reflection light intensity.
  float3          P;                                                                                // Scene position.
  float3          N;                                                                                // Scene normal direction.
  
  // INITIALIZING RAY MARCHING:
  V[0] = V_mat[0];                                                                                  // Getting view matrix...
  V[1] = V_mat[1];                                                                                  // Getting view matrix...
  V[2] = V_mat[2];                                                                                  // Getting view matrix...
  V[3] = V_mat[3];                                                                                  // Getting view matrix...

  camera.pos = (float3)(0.0f, 0.0f, 0.0f);                                                          // Setting initial camera position...
  camera.fov = FOV;                                                                                 // Setting camera FOV...

  light.pos = light_position[0].xyz;                                                                // Setting light position...
  light.col = light_color[0];                                                                       // Setting light color [r, g, b, ambient]...

  M.amb = material[0];                                                                              // Setting material ambient color...
  M.dif = material[1];                                                                              // Setting material diffusion color...
  M.ref = material[2];                                                                              // Setting material reflection color...
  M.shn = material[3];                                                                              // Setting material shininess...

  camera.pos = mul(V, (float4)(camera.pos, 1.0f)).xyz;                                              // Applying arcball to camera position...
  ray = (x, y, -1.0f/tan(camera.fov*PI/360.0f));                                                    // Computing ray intersection on canvas...
  ray = mul(V, (float4)(ray, 1.0f)).xyz;                                                            // Applying arcball to ray direction...
  ray = normalize(ray);                                                                             // Normalizing ray direction...
  
  // COMPUTING RAY MARCHING:
  d = raymarch(camera.pos, ray);                                                            // Computing scene distance...
  P = camera.pos + d*ray;                                                                   // Computing scene position...
  N = normal(P);                                                                                    // Computing scene normal...
  
  // COMPUTING LIGHTNING:
  view = normalize(camera.pos - P);                                                                 // Computing view direction...
  incident = normalize(light.pos - P);                                                              // Computing incident light direction...
  reflected = reflect(-incident, N);                                                                // Computing reflected light direction...
  halfway = normalize(incident + view);                                                             // Coputing halfway vector (Blinn-Phong)...
  ambient = light.col.a;                                                                            // Computing light ambient color...
  reflection = pow(max(dot(N, halfway), 0.0f), length(M.shn));                                      // Computing light reflection intensity
  diffusion = clamp(dot(N, incident), 0.0f, 1.0f)*shadow(P + N*2.0f*EPSILON, incident, 10.0f);      // Computing light diffusion intensity... 
  
  fragment.amb = (float4)(ambient*light.col.xyz*M.amb.xyz, 0.0f);                                                           // Computing light ambient color...
  fragment.dif = (float4)(diffusion*light.col.xyz*M.dif.xyz, 0.0f);                                                         // Computing light diffused color...
  fragment.ref = (float4)(reflection*light.col.xyz*M.ref.xyz, 0.0f);                                                        // Computing light reflected color...

  fragment_color[i + (uint)W*j] = fragment.amb + fragment.dif + fragment.ref;                       // Setting output color...
  fragment_color[i + (uint)W*j].a = 1.0f;
}
