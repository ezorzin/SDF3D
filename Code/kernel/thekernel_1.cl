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
  float3 pos;                                                                                        // Camera position.
  float  fov;                                                                                        // Camera field of view (degrees).
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
  
  return normalize((float3)(nx, ny, nz));                                                           // Returning normal vector...
}

float3 reflect(float3 I, float3 N)
{
  float3 reflect;

  reflect = I - 2.0f*dot(N, I)*N;

  return reflect;
}

__kernel void thekernel(__global float4*    fragment_color,                                         // Fragment color.
                        __global float4*    V,                                                      // View matrix [4x4].
                        __global float4*    canvas_param,                                           // Canvas [W, H, AR, FOV].
                        __global float4*    camera_param,                                           // Camera [x, y, z, 1.0f].
                        __global float4*    light_param,                                            // Light  [x, y, z, 1.0f].
                        __global float4*    material_param                                          // Material.
                        )
{ 
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////// INDICES //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  uint         i = get_global_id(0);                                                                // Global index i [#].
  uint         j = get_global_id(1);                                                                // Global index j [#].

  float        W = canvas_param[0].x;
  float        H = canvas_param[0].y;
  float        AR = canvas_param[0].z;
  float        FOV = canvas_param[0].w;

  float        x = 2.0f*(i/W) - 1.0f;
  float        y = 2.0f*(j/H) - 1.0f;

  float4       V_0 = V[0];
  float4       V_1 = V[1];
  float4       V_2 = V[2];
  float4       V_3 = V[3];

  float4       camera  = camera_param[0];
  float        camera_x;
  float        camera_y;
  float        camera_z;
  float        camera_w;

  float4       camera_pos = camera_param[0];                                                        // Camera position [x, y, z, w].
  float4       light_pos = light_param[0];                                                          // Light position [x, y, z, w].
  float4       light_col = light_param[1];                                                          // Light color [r, g, b, ambient].
  float4       M_amb = material_param[0];
  float4       M_dif = material_param[1];
  float4       M_ref = material_param[2];
  float4       M_shn = material_param[3];
  

  float4      ray;                                                                                  // Marching ray.
  float       ray_x;
  float       ray_y;
  float       ray_z;
  float       ray_w;

  float       d;                                                                                    // Marching distance.
  float3      view;                                                                                 // View direction.
  float3      incident;                                                                             // Incident light direction.
  float3      reflected;                                                                            // Reflected light direction.
  float3      halfway;                                                                              // Halfway light direction (Blinn-Phong).
  float4      ambient;                                                                              // Ambient light color.
  float4      diffusion;                                                                            // Diffusion light color.
  float4      reflection;                                                                           // Reflection light color.
  float3      P;                                                                                    // Scene position.
  float3      N;                                                                                    // Scene normal direction.
  float       light_dif;
  float       light_ref;

  // INITIALIZING RAY MARCHING:
  camera_x = dot(V_0, camera);                                                                      // Applying arcball to camera position...
  camera_y = dot(V_1, camera);                                                                      // Applying arcball to camera position...
  camera_z = dot(V_2, camera);                                                                      // Applying arcball to camera position...
  camera_w = dot(V_3, camera);                                                                      // Applying arcball to camera position...
  camera.x = camera_x;
  camera.y = camera_y;
  camera.z = camera_z;
  camera.w = camera_w;                                  

  //printf("%f\n", camera.w);

  ray = (float4)(x*AR, y, -2.0f/tan(FOV*PI/360.0f), 1.0f);                                          // Computing position on canvas (ray intersection on quad)...
  ray_x = dot(V_0, ray);                                                                            // Applying arcball to ray direction...
  ray_y = dot(V_1, ray);                                                                            // Applying arcball to ray direction...
  ray_z = dot(V_2, ray);                                                                            // Applying arcball to ray direction...
  ray_w = dot(V_3, ray);                                                                            // Applying arcball to ray direction...
  ray = normalize((float4)(ray_x, ray_y, ray_z, 0.0f));                                            
  
  // COMPUTING RAY MARCHING:
  d = raymarch(camera_pos.xyz, ray.xyz);                                                            // Computing scene distance...
  P = camera_pos.xyz + d*ray.xyz;                                                                   // Computing scene position...
  N = normal(P);                                                                                    // Computing scene normal...
  
  // COMPUTING LIGHTNING:
  view = normalize(camera_pos.xyz - P);                                                             // Computing view direction...
  incident = normalize(light_pos.xyz - P);                                                          // Computing incident light direction...
  reflected = reflect(-incident, N);                                                                // Computing reflected light direction...
  halfway = normalize(incident + view);                                                             // Coputing halfway vector (Blinn-Phong)...
  light_ref = pow(max(dot(N, halfway), 0.0f), length(M_shn));                                       // Computing light reflection intensity
  light_dif = clamp(dot(N, incident), 0.0f, 1.0f)*shadow(P + N*2.0f*EPSILON, incident, 10.0f);      // Computing light diffusion intensity... 
  ambient = light_col.a*M_amb;                                                                      // Computing light ambient color...
  diffusion = light_dif*M_dif;                                                                      // Computing light diffused color...
  reflection = light_ref*M_ref;                                                                     // Computing light reflected color...

  fragment_color[i + (uint)W*j] = ambient + diffusion + reflection;                                 // Setting output color...
  fragment_color[i + (uint)W*j].a = 1.0f;
}
