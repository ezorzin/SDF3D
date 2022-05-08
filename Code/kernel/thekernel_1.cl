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
  float3  col;                                                                                      // Light color.
  float   amb;                                                                                      // Light ambient intensity.
  float   k;                                                                                        // Light sharpness.
};

// BALL:
struct Ball
{
  float3  pos;                                                                                      // Ball position.
  float   rad;                                                                                      // Ball radius.
};

// MATERIAL:
struct Material
{
  float3  amb;                                                                                      // Material ambient color.
  float3  dif;                                                                                      // Material diffused color.
  float3  ref;                                                                                      // Material reflected color.
  float   t;                                                                                        // Material transparency.
  float   n;                                                                                        // Material index of refraction.
  float   s;                                                                                        // Material shininess.
};

// FRAGMENT:
struct Fragment
{
  float3  amb;
  float3  dif;
  float3  ref;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// OBJECT'S SIGNED DISTANCE FIELDS IMPLICIT FUNCTIONS /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
float sphereSDF(float3 ray, struct Ball ball)
{
  float sdf = length(ray - ball.pos) - ball.rad;                                                    // Signed distance field.
 
  return sdf;                                                                                       // Returning signed distance field...
}

float planeSDF(float3 ray)
{
  float sdf = ray.y;                                                                                // Signed distance field.
 
  return sdf;                                                                                       // Returning signed distance field...
}

float sceneSDF(float3 ray, struct Ball ball[], int n)
{
  float sdf = INF;                                                                                  // Signed distance field.
  int i;

  sdf = min(sdf, planeSDF(ray));                                                                    // Computing signed distance field...
  
  for (i = 0; i < n; i++)
  {
    sdf = min(sdf, sphereSDF(ray, ball[i]));                                                             // Computing signed distance field...
  }

  return sdf;                                                                                       // Returning signed distance field...
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// RAY MARCHING ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
float raymarch(float3 position, float3 direction, struct Ball ball[], int n)
{
  float   distance = 0.0f;                                                                          // Marching distance.
  float   sdf;                                                                                      // Signed distance field.
  int     i;                                                                                        // Step index.
  float3  ray;                                                                                      // Marching ray.

  for (i = 0; i < MAX_STEPS; i++)
  {
    ray = position + distance*direction;                                                            // Computing marching ray...    
    sdf = sceneSDF(ray, ball, n);                                                                            // Computing scene distance field...
    distance += sdf;                                                                                // Updating marching distance...

    if(distance > MAX_DISTANCE || sdf < EPSILON) break;                                             // Checking numeric precision constraints...
  }

  return distance;                                                                                  // Returning final marching distance...
}

float shadow(float3 position, float3 direction, struct Ball ball[], int n, struct Light light)
{
  float   distance = 0.0f;                                                                          // Marching distance.
  float   sdf = INF;                                                                                // Previous step signed distance field.
  float   sdf_new;                                                                                  // Current step signed distance field.
  int     i;                                                                                        // Step index.
  float3  ray;                                                                                      // Marching ray.
  float   shadow = 1.0f;                                                                            // Shadow intensity.
  float   intersection;                                                                             // Previous to current SDF on-ray intersection point.
  float   d_est;                                                                                    // Estimated closest distance.

  for (i = 0; i < MAX_STEPS; i++)
  {
    ray = position + distance*direction;                                                            // Computing marching ray...   
    sdf_new = sceneSDF(ray, ball, n);                                                                  // Computing scene distance field... 
    intersection = (i == 0) ? 0.0f : sdf_new*sdf_new/(2.0f*sdf);                                    // Computing on-ray intersection point...
    d_est = sqrt(sdf_new*sdf_new - intersection*intersection);                                      // Computing estimated closest distance...
    shadow = min(shadow, light.k*d_est/max(0.0f, distance - intersection));                         // Computing shadow intensity...
    sdf = sdf_new;                                                                                  // Backing up signed distance field...
    distance += sdf_new;                                                                            // Updating marching distance...
    
    if(distance > MAX_DISTANCE || shadow < EPSILON) break;                                          // Checking numeric precision constraints...
  }

  shadow = clamp(shadow, 0.0f, 1.0f);                                                               // Clamping shadow intensity...

  return shadow;                                                                                    // Returning final shadow intensity...
}

float3 normal(float3 position, struct Ball ball[], int n)
{
  float nx;                                                                                         // Scene SDF's gradient x-component.
  float ny;                                                                                         // Scene SDF's gradient y-component.
  float nz;                                                                                         // Scene SDF's gradient z-component.
  float sdf_L;                                                                                      // Scene SDF (left limit).
  float sdf_R;                                                                                      // Scene SDF (right limit).

  sdf_L = sceneSDF(position - DX, ball, n);                                                                  // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DX, ball, n);                                                                  // Computing signed distance field (left limit)...
  nx = sdf_R - sdf_L;                                                                               // Computing gradient (x-component, centered derivative)...

  sdf_L = sceneSDF(position - DY, ball, n);                                                                  // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DY, ball, n);                                                                  // Computing signed distance field (left limit)...
  ny = sdf_R - sdf_L;                                                                               // Computing gradient (x-component, centered derivative)...

  sdf_L = sceneSDF(position - DZ, ball, n);                                                                  // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DZ, ball, n);                                                                  // Computing signed distance field (left limit)...
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

__kernel void thekernel(__global float4*    fragment_color,                                         // Fragment color.
                        __global float4*    center,
                        __global float4*    view_matrix,                                            // View matrix [4x4].
                        __global float4*    canvas,                                                 // Canvas [W, H, AR, FOV].
                        __global float4*    light_position,                                         // Light position [x, y, z, k].
                        __global float4*    light_color,                                            // Light color [r, g, b, ambient].
                        __global float4*    ball_position,                                          // Ball position [x, y, z, radius].
                        __global float4*    material_ambient,                                       // Material ambient color [r, g, b, transparency].
                        __global float4*    material_diffusion,                                     // Material diffusion color [r, g, b, n_index].
                        __global float4*    material_reflection                                     // Material reflection color [r, g, b, shininess].
                        )
{ 
  uint            i = get_global_id(0);                                                             // Global index i [#].
  uint            j = get_global_id(1);                                                             // Global index j [#].
  uint            k;

  float4          V[4];                                                                             // View matrix.

  float           W = canvas[0].x;                                                                  // Window width [px].
  float           H = canvas[0].y;                                                                  // Window height [px].
  float           AR = canvas[0].z;                                                                 // Window aspect ratio [].
  float           FOV = canvas[0].w;                                                                // FOV [degrees].
  float           x = 2.0f*(i/W) - 1.0f;                                                            // Canvas x-coordinate [-1.0f...+1.0f].
  float           y = 2.0f*(j/H) - 1.0f;                                                            // Canvas y-coordinate [-1.0f...+1.0f].

  struct Camera   camera;                                                                           // Camera.
  struct Light    light;                                                                            // Light.
  struct Ball     ball[3];                                                                             // Ball.
  struct Material material;                                                                         // Material.
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
  V[0] = view_matrix[0];                                                                            // Getting view matrix...
  V[1] = view_matrix[1];                                                                            // Getting view matrix...
  V[2] = view_matrix[2];                                                                            // Getting view matrix...
  V[3] = view_matrix[3];                                                                            // Getting view matrix...

  camera.pos = (float3)(0.0f, 0.2f, 2.0f);                                                          // Setting initial camera position...
  camera.fov = FOV;                                                                                 // Setting camera FOV...

  light.pos = light_position[0].xyz;                                                                // Setting light position...
  light.col = light_color[0].xyz;                                                                   // Setting light color [r, g, b]...
  light.amb = light_color[0].w;                                                                     // Setting light ambient intensity [ambient]...
  light.k   = light_position[0].w;                                                                  // Setting light sharpness [k]...

  for (k = 0; k < 3; k++)
  {
    ball[k].pos = ball_position[k].xyz;                                                                  // Setting ball position [x, y, z]...
    ball[k].rad = ball_position[k].w;                                                                    // Setting ball radius [radius]...
  }

  material.amb = material_ambient[0].xyz;                                                           // Setting material ambient color...
  material.dif = material_diffusion[0].xyz;                                                         // Setting material diffusion color...
  material.ref = material_reflection[0].xyz;                                                        // Setting material reflection color...
  material.t   = material_ambient[0].w;                                                             // Setting materila transparency...
  material.n   = material_diffusion[0].w;                                                           // Setting material index of refraction...
  material.s   = material_reflection[0].w;                                                          // Setting material shininess...

  camera.pos = mul(V, (float4)(camera.pos, 1.0f)).xyz;                                              // Applying arcball to camera position...
  ray = (float3)(x*AR, y, -1.0f/tan(camera.fov*PI/360.0f));                                                    // Computing ray intersection on canvas...
  ray = mul(V, (float4)(ray, 1.0f)).xyz;                                                            // Applying arcball to ray direction...
  ray = normalize(ray);                                                                             // Normalizing ray direction...
  
  
  // COMPUTING RAY MARCHING:
  d = raymarch(camera.pos, ray, ball, 3);                                                                    // Computing scene distance...
  P = camera.pos + d*ray;                                                                           // Computing scene position...
  N = normal(P, ball, 3);                                                                                    // Computing scene normal...
  
  // COMPUTING LIGHTNING:
  view = normalize(camera.pos - P);                                                                 // Computing view direction...
  incident = normalize(light.pos - P);                                                              // Computing incident light direction...
  reflected = reflect(-incident, N);                                                                // Computing reflected light direction...
  halfway = normalize(incident + view);                                                             // Coputing halfway vector (Blinn-Phong)...
  ambient = light.amb;                                                                              // Computing light ambient color...
  reflection = pow(max(dot(N, halfway), 0.0f), material.s);                                         // Computing light reflection intensity
  diffusion = clamp(dot(N, incident), 0.0f, 1.0f)*shadow(P + N*2.0f*EPSILON, incident, ball, 3, light);// Computing light diffusion intensity... 
  
  fragment.amb = ambient*light.col*material.amb;                                                    // Computing light ambient color...
  fragment.dif = diffusion*light.col*material.dif;                                                  // Computing light diffused color...
  fragment.ref = reflection*light.col*material.ref;                                                 // Computing light reflected color...

  fragment_color[i + (uint)W*j] = (float4)(fragment.amb + fragment.dif + fragment.ref, 1.0f);       // Setting output color...
}
