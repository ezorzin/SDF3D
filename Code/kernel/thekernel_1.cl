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

// FRAGMENT:
struct Fragment
{
  float3  amb;
  float3  dif;
  float3  ref;
};

enum ObjectType
{
  SCENE = 0,
  PLANE = 1,
  SPHERE = 2
}; 

// Object:
struct Object
{
  int     type;                                                                                     // Object type.
  float16 T;                                                                                        // Object transformation matrix.
  float4  par;                                                                                      // Object parameters.
  float4  amb;                                                                                      // Object (ambient color, transparency).
  float4  dif;                                                                                      // Object (diffused color, shininess).
  float4  ref;                                                                                      // Object (reflected color, index of refraction).
  float   sdf;                                                                                      // Signed distance field.
};

float4 mul(float16 M, float4 V)
{
  float4 A;

  A.x = M.s0*V.x + M.s1*V.y + M.s2*V.z + M.s3*V.w;
  A.y = M.s4*V.x + M.s5*V.y + M.s6*V.z + M.s7*V.w;
  A.z = M.s8*V.x + M.s9*V.y + M.sA*V.z + M.sB*V.w;
  A.w = M.sC*V.x + M.sD*V.y + M.sE*V.z + M.sF*V.w;

  return A;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////// OBJECT'S SIGNED DISTANCE FIELDS IMPLICIT FUNCTIONS /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
float objectSDF (int type, float16 T, float4 parameter, float3 ray)
{
  float sdf;                                                                                        // Signed distance field.

  ray = mul(T, (float4)(ray, 1.0f)).xyz;                                                            // Applying transformation matrix...

  switch (type)
  {
    case SCENE:
      sdf = INF;

    case PLANE:
      sdf = dot(ray, (float3)(0.0f, 0.0f, 1.0f));                                                   // Computing sdf...
      break;

    case SPHERE:
      // parameter.x = radius.
      sdf = length(ray) - parameter.x;                                                              // Computing sdf...
      break;
  }
 
  return sdf;                                                                                       // Returning sdf...
}

struct Object unionSDF(struct Object a, struct Object b)
{
  return a.sdf < b.sdf? a : b;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// RAY MARCHING ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
struct Object raymarch(float3 position, float3 direction, struct Object object)
{
  float             distance = 0.0f;                                                                // Marching distance.
  struct Object     scene;                                                                          // Scene.
  int               i;                                                                              // Step index.
  float3            ray;                                                                            // Marching ray.

  for (i = 0; i < MAX_STEPS; i++)
  {
    ray = position + distance*direction;                                                            // Computing marching ray...    
    distance += object.sdf;                                                                         // Updating marching distance...

    if(distance > MAX_DISTANCE || object.sdf < EPSILON)                                             // Checking numeric precision constraints...
    {
      scene = object;
      break;
    }
  }

  scene.sdf = distance;                                                                             // Setting scene sdf...

  return scene;                                                                                     // Returning scene (color, sdf)...
}

float shadow(float3 position, float3 direction, struct Object object, struct Light light)
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
    sdf_new = object.sdf;                                                                           // Setting distance field... 
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

  sdf_L = sceneSDF(position - DX, ball, n).w;                                                       // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DX, ball, n).w;                                                       // Computing signed distance field (left limit)...
  nx = sdf_R - sdf_L;                                                                               // Computing gradient (x-component, centered derivative)...

  sdf_L = sceneSDF(position - DY, ball, n).w;                                                       // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DY, ball, n).w;                                                       // Computing signed distance field (left limit)...
  ny = sdf_R - sdf_L;                                                                               // Computing gradient (x-component, centered derivative)...

  sdf_L = sceneSDF(position - DZ, ball, n).w;                                                       // Computing signed distance field (left limit)...
  sdf_R = sceneSDF(position + DZ, ball, n).w;                                                       // Computing signed distance field (left limit)...
  nz = sdf_R - sdf_L;                                                                               // Computing gradient (x-component, centered derivative)...   
  
  return normalize((float3)(nx, ny, nz));                                                           // Returning normal vector...
}

float3 reflect(float3 I, float3 N)
{
  float3 reflect;

  reflect = I - 2.0f*dot(N, I)*N;

  return reflect;
}

float3 refract(float3 I, float3 N, float n1, float n2)
{
  float n;
  float k;
  float3 refract;

  n = n1/n2;
  k = 1.0f - n*n*(1.0f - dot(N, I)*dot(N, I));
  
  if (k < 0.0f)
  {
    refract = (float3)(0.0f, 0.0f, 0.0f);
  }
  else
  {
    refract = n*I - (n*dot(N, I) + sqrt(k))*N;
  }
        
  return refract;  
}

__kernel void thekernel(__global float4*    fragment_color,                                         // Fragment color.
                        __global float4*    center,
                        __global float16*   view_matrix,                                            // View matrix [4x4].
                        __global float4*    canvas,                                                 // Canvas [W, H, AR, FOV].
                        __global float4*    light_position,                                         // Light position [x, y, z, k].
                        __global float4*    light_color,                                            // Light color [r, g, b, ambient].
                        __global int*       object_type,                                            // Object type.
                        __global float16*   object_transformation,                                  // Object transformation matrix.                     
                        __global float4*    object_ambient,                                         // Object ambient color [r, g, b, transparency].
                        __global float4*    object_diffusion,                                       // Object diffusion color [r, g, b, n_index].
                        __global float4*    object_reflection,                                      // Object reflection color [r, g, b, shininess].
                        __global float4*    object_parameter                                        // Object parameters.
                        )
{ 
  uint            i = get_global_id(0);                                                             // Global index i [#].
  uint            j = get_global_id(1);                                                             // Global index j [#].
  uint            k;
  uint            offset = 0;

  uint num_objects; // EZOR: to be defined in kernel arguments.

  float16         M;                                                                                // View matrix.
  float16         T;                                                                                // Tranformation matrix.

  float           W = canvas[0].x;                                                                  // Window width [px].
  float           H = canvas[0].y;                                                                  // Window height [px].
  float           AR = canvas[0].z;                                                                 // Window aspect ratio [].
  float           FOV = canvas[0].w;                                                                // FOV [degrees].
  float           x = 2.0f*(i/W) - 1.0f;                                                            // Canvas x-coordinate [-1.0f...+1.0f].
  float           y = 2.0f*(j/H) - 1.0f;                                                            // Canvas y-coordinate [-1.0f...+1.0f].

  struct Camera   camera;                                                                           // Camera.
  struct Light    light;                                                                            // Light.
  struct Object   object;                                                                           // Object.
  struct Object   scene;                                                                            // Scene.
  struct Fragment fragment;                                                                         // Fragment.

  float3          ray;                                                                              // Marching ray.
  float           distance;                                                                         // Marching distance.
  float3          view;                                                                             // View direction.
  float3          incident;                                                                         // Incident light direction.
  float3          reflected;                                                                        // Reflected light direction.
  float3          refracted;                                                                        // Refracted light direction.
  float3          halfway;                                                                          // Halfway light direction (Blinn-Phong).
  float           ambient;                                                                          // Ambient light intensity.
  float           diffusion;                                                                        // Diffusion light intensity.
  float           reflection;                                                                       // Reflection light intensity.
  float3          P;                                                                                // Scene position.
  float3          N;                                                                                // Scene normal direction.

  // INITIALIZING RAY MARCHING:
  M = view_matrix[0];                                                                               // Getting view matrix...

  camera.pos = (float3)(0.0f, 0.2f, 2.0f);                                                          // Setting initial camera position...
  camera.pos = mul(M, (float4)(camera.pos, 1.0f)).xyz;                                              // Applying arcball to camera position...
  camera.fov = FOV;                                                                                 // Setting camera FOV...

  light.pos = light_position[0].xyz;                                                                // Setting light position...
  light.col = light_color[0].xyz;                                                                   // Setting light color [r, g, b]...
  light.amb = light_color[0].w;                                                                     // Setting light ambient intensity [ambient]...
  light.k   = light_position[0].w;                                                                  // Setting light sharpness [k]...

  ray = (float3)(x*AR, y, -1.0f/tan(camera.fov*PI/360.0f));                                         // Computing ray intersection on canvas...
  ray = mul(M, (float4)(ray, 1.0f)).xyz;                                                            // Applying arcball to ray direction...
  ray = normalize(ray);                                                                             // Normalizing ray direction...

  // COMPUTING RAY MARCHING:
  distance = 0.0f;                                                                                  // Marching distance.
  
  for (i = 0; i < MAX_STEPS; i++)
  {
    P = camera.pos + distance*ray;                                                                  // Computing marching ray...    

    scene.type = 0;                                                                                   // Setting object type...
    scene.T = 0.0f;                                                                                   // Setting object transformation matrix...
    scene.par = 0.0f;                                                                                 // Setting object parameters...
    scene.amb = 0.0f;                                                                                 // Setting object ambient color...
    scene.dif = 0.0f;                                                                                 // Setting object diffusion color...
    scene.ref = 0.0f;                                                                                 // Setting object reflection color...
    scene.sdf = INF;

    for (k = 0; k < num_objects; k++)
    {
      object.type = object_type[k];                                                                 // Getting object type...
      object.T = object_transformation[k];                                                          // Getting object transformation matrix...
      object.par = object_parameter[k];                                                             // Getting object parameters...
      //object.amb = object_ambient[k];                                                               // Getting object ambient color...
      //object.dif = object_diffusion[k];                                                             // Getting object diffusion color...
      //object.ref = object_reflection[k];                                                            // Getting object reflection color...
      sdf =  objectSDF (object_type[k], object_transformation[k], object_parameter[k], P);          // Computing object sdf...
                                                                  
      scene = unionSDF(scene, object);                                                              // Assembling scene...
    }

    distance += scene.sdf;                                                                          // Updating marching distance...

    if(distance > MAX_DISTANCE || scene.sdf < EPSILON)                                              // Checking numeric precision constraints...
    {
      P = camera.pos + distance*ray;                                                                // Computing scene position...
      scene.sdf = distance;                                                                         // Setting scene sdf...
      break;
    }
  }

  

  N = normal(P, ball, 3);                                                                           // Computing scene normal...

  




  
  
  
  
  
  col_d = raymarch(camera.pos, ray, ball, 3);                                                       // Computing scene distance...
 
  
  // COMPUTING LIGHTNING:
  view = normalize(camera.pos - P);                                                                 // Computing view direction...
  incident = normalize(light.pos - P);                                                              // Computing incident light direction...
  
  reflected = reflect(-incident, N);                                                                // Computing reflected light direction...
  halfway = normalize(incident + view);                                                             // Coputing halfway vector (Blinn-Phong)...
  ambient = light.amb;                                                                              // Computing light ambient color...
  reflection = pow(max(dot(N, halfway), 0.0f), material.s);                                         // Computing light reflection intensity
  diffusion = clamp(dot(N, incident), 0.0f, 1.0f)*shadow(P + N*2.0f*EPSILON, incident, ball, 3, light);// Computing light diffusion intensity... 
  
  /*
  refracted = refract(incident, N, 1.00f, 1.33f);
  d = raymarch(P, -refracted, ball, 3).w;                                                                    // Computing scene distance...
  P = P + d*ray;                                                                           // Computing scene position...
  N = normal(P, ball, 3);                                                                           // Computing scene normal...
  */

  fragment.amb = ambient*light.col*material.amb;                                                    // Computing light ambient color...
  fragment.dif = diffusion*light.col*col_d.xyz;                                                  // Computing light diffused color...
  fragment.ref = reflection*light.col*material.ref;                                                 // Computing light reflected color...

  fragment_color[i + (uint)W*j] = (float4)(fragment.amb + fragment.dif + fragment.ref, 1.0f);       // Setting output color...
}
