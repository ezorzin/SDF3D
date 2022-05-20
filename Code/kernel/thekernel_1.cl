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

#define OOO  0.0f, 0.0f, 0.0f
#define III  1.0f, 1.0f, 1.0f
#define IOO  1.0f, 0.0f, 0.0f
#define OIO  0.0f, 1.0f, 0.0f
#define OOI  0.0f, 0.0f, 1.0f
#define ZERO3 (float3)(OOO)
#define ONE3  (float)(III)

#define TET_A (float3)(+1.0f, -1.0f, -1.0f)
#define TET_B (float3)(-1.0f, -1.0f, +1.0f)
#define TET_C (float3)(-1.0f, +1.0f, -1.0f)
#define TET_D (float3)(+1.0f, +1.0f, +1.0f)

#define OOOO 0.0f, 0.0f, 0.0f, 0.0f
#define IIII 1.0f, 1.0f, 1.0f, 1.0f
#define IOOO 1.0f, 0.0f, 0.0f, 0.0f
#define OIOO 0.0f, 1.0f, 0.0f, 0.0f
#define OOIO 0.0f, 0.0f, 1.0f, 0.0f
#define OOOI 0.0f, 0.0f, 0.0f, 1.0f
#define ZERO4 (float4)(OOOO)
#define I4x4 (float16)(IOOO, OIOO, OOIO, OOOI)

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////// DATA STRUCTURES //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
struct Camera
{
  float3 pos;                                                                                       // Camera position.
  float  fov;                                                                                       // Camera field of view (degrees).
};

struct Light
{
  float3  pos;                                                                                      // Light position.
  float3  col;                                                                                      // Light color.
  float   amb;                                                                                      // Light ambient intensity.
  float   k;                                                                                        // Light sharpness.
};

struct Fragment
{
  float3  amb;                                                                                      // Fragment ambient light.
  float3  dif;                                                                                      // Fragment diffused light.
  float3  ref;                                                                                      // Fragment reflected light.
};

struct Object
{
  int     type;                                                                                     // Object type.
  float16 T;                                                                                        // Object transformation matrix.
  float4  par;                                                                                      // Object parameters.
  float4  amb;                                                                                      // Object (ambient color, transparency).
  float4  dif;                                                                                      // Object (diffused color, shininess).
  float4  ref;                                                                                      // Object (reflected color, index of refraction).
  float   shd;                                                                                      // Object shadow.
  float   sdf;                                                                                      // Signed distance field.
  float3  P;                                                                                        // Current point.
  float3  N;                                                                                        // Current normal.
};

enum ObjectType
{
  SCENE = 0,
  PLANE = 1,
  SPHERE = 2
}; 

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// LINEAR ALGEBRA //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////// SIGNED DISTANCE FIELD ///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
float3 displaceSDF(struct Object object, float3 point)
{
  point = mul(object.T, (float4)(point, 1.0f)).xyz;                                                            // Applying transformation matrix...

  return point;
}

float objectSDF (struct Object object, float3 point)
{
  float sdf;

  point = displaceSDF(object, point);

  switch (object.type)
  {
    case SCENE:
      sdf = INF;

    case PLANE:
      sdf = dot(point, (float3)(0.0f, 1.0f, 0.0f));                                                   // Computing sdf...
      sdf = point.y;
      break;

    case SPHERE:
      // parameter.x = radius.
      sdf = length(point) - object.par.x;                                                       // Computing sdf...
      break;
  }
 
  return sdf;                                                                                       // Returning object...
}

struct Object unionSDF(struct Object a, struct Object b)
{
  return a.sdf < b.sdf? a : b;
}

float3 normalSDF(struct Object object, float3 point, float h)
{
  float3 N;
  int i;

  N = ZERO3;
  
  N += objectSDF(object, point + TET_A*h)*TET_A;
  N += objectSDF(object, point + TET_B*h)*TET_B;
  N += objectSDF(object, point + TET_C*h)*TET_C;
  N += objectSDF(object, point + TET_D*h)*TET_D;

  return normalize(N);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// RAY MARCHING ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
struct Object raymarchSDF(struct Object object, float3 origin, float3 direction, float h)
{
  float             distance = 0.0f;                                                                // Marching distance.
  int               i;                                                                              // Step index.
  int               k;                                                                              // Object index.
  float3            P;                                                                              // Marching ray.
  float3            P_mod;                                                                          // Modified marching ray.

  for (i = 0; i < MAX_STEPS; i++)
  {
    P = origin + distance*direction;                                                                // Computing marching ray...
    object.sdf = objectSDF(object, P);                                                              // Computing object sdf...
    distance += object.sdf;                                                                         // Updating marching distance...
    
    if(distance > MAX_DISTANCE || object.sdf < EPSILON) break;                                      // Checking numeric precision constraints...
  }

  if(distance > MAX_DISTANCE)
  {
    object.amb = ZERO4;
    object.dif = ZERO4;
    object.ref = ZERO4;
  }

  object.sdf = distance;                                                                            // Computing object distance...
  P = origin + distance*direction;                                                                  // Computing marching ray...
  object.P = P;                                                                                     // Computing object ray position...                                            
  object.N = normalSDF(object, P, h);                                                               // Computing object normal...
  
  return object;                                                                                    // Returning scene...
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// OPTICS ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
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
/*
                        __global float4*    object_parameter                                        // Object parameters.
                        __global float4*    object_ambient,                                         // Object ambient color [r, g, b, transparency].
                        __global float4*    object_diffusion,                                       // Object diffusion color [r, g, b, n_index].
                        __global float4*    object_reflection,                                      // Object reflection color [r, g, b, shininess].
                        
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// KERNEL ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
__kernel void thekernel(__global float4*    fragment_color,                                         // Fragment color.
                        __global float4*    center,
                        __global float16*   view_matrix,                                            // View matrix [4x4].
                        __global float4*    canvas,                                                 // Canvas [W, H, AR, FOV].
                        __global float4*    light_position,                                         // Light position [x, y, z, k].
                        __global float4*    light_color,                                            // Light color [r, g, b, ambient].
                        __global int*       object_number,                                          // Object number.
                        __global int*       object_type,                                            // Object type.
                        __global float16*   T,                                                      // Object transformation matrix [4x4].                     
                        __global float16*   M                                                       // Object material matrix [4x4].
                        )
{ 
  uint            i = get_global_id(0);                                                             // Global index i [#].
  uint            j = get_global_id(1);                                                             // Global index j [#].
  uint            k;                                                                                // Object index [#].
  int             n = object_number[0];                                                             // Number of object [#].
  float16         V;                                                                                // View matrix [4x4].
  float           W = canvas[0].x;                                                                  // Window width [px].
  float           H = canvas[0].y;                                                                  // Window height [px].
  float           AR = canvas[0].z;                                                                 // Window aspect ratio [].
  float           FOV = canvas[0].w;                                                                // FOV [degrees].
  float           x = 2.0f*(i/W) - 1.0f;                                                            // Canvas x-coordinate [-1.0f...+1.0f].
  float           y = 2.0f*(j/H) - 1.0f;                                                            // Canvas y-coordinate [-1.0f...+1.0f].
  
  struct Camera   camera;                                                                           // Camera.
  struct Light    light;                                                                            // Light.
  struct Object   object;                                                                            // Scene.
  struct Object   object2;                                                                            // Scene.
  struct Object   scene;                                                                            // Scene.
  struct Object   scene2;                                                                            // Scene.
  struct Fragment fragment;                                                                         // Fragment.

  float3          ray;                                                                              // Marching ray.
  float           distance;                                                                         // Marching distance.
  float3          view;                                                                             // View direction.
  float3          incident;                                                                         // Incident light direction.
  float3          incident2;                                                                         // Incident light direction.
  float3          reflected;                                                                        // Reflected light direction.
  float3          refracted;                                                                        // Refracted light direction.
  float3          halfway;                                                                          // Halfway light direction (Blinn-Phong).
  float           ambient;                                                                          // Ambient light intensity.
  float           diffusion;                                                                        // Diffusion light intensity.
  float           reflection;                                                                       // Reflection light intensity.
  float           shadow;                                                                           // Shadow intensity.
  float           shadow_old;
  float3          P;                                                                                // Scene position.
  float3          N;                                                                                // Scene normal direction.
  float           h;                                                                                // Normal precision (for antialiasing).
  bool first = true;

  float d[6];

  // INITIALIZING RAY MARCHING:
  V = view_matrix[0];                                                                               // Getting view matrix...

  camera.pos = (float3)(0.0f, 0.2f, 2.0f);                                                          // Setting initial camera position...
  camera.pos = mul(V, (float4)(camera.pos, 1.0f)).xyz;                                              // Applying arcball to camera position...
  camera.fov = FOV;                                                                                 // Setting camera FOV...

  light.pos = light_position[0].xyz;                                                                // Setting light position...
  light.k   = light_position[0].w;                                                                  // Setting light sharpness [k]...
  light.col = light_color[0].xyz;                                                                   // Setting light color [r, g, b]...
  light.amb = light_color[0].w;                                                                     // Setting light ambient intensity [ambient]...

  ray = (float3)(x*AR, y, -1.0f/tan(camera.fov*PI/360.0f));                                         // Computing ray intersection on canvas...
  ray = mul(V, (float4)(ray, 1.0f)).xyz;                                                            // Applying arcball to ray direction...
  ray = normalize(ray);                                                                             // Normalizing ray direction...

  h = EPSILON;                                                                                       // EZOR: to be better defined...

  n = 2;

  scene.sdf = INF;
  scene.shd = 1.0f;

  // COMPUTING RAY MARCHING:
  for(k = 0; k < n; k++)
  {
    object.type = object_type[k];                                                                   // Getting object type...
    object.T = T[k];                                                                                // Getting object transformation matrix...
    object.par = M[k].s0123;                                                                        // Getting object parameters...
    object.amb = M[k].s4567;                                                                        // Getting object ambient color...
    object.dif = M[k].s89AB;                                                                        // Getting object diffusion color...
    object.ref = M[k].sCDEF;                                                                        // Getting object reflection color...
    
    object = raymarchSDF(object, camera.pos, ray, h);                                               // Computing object raymarching...

    object2 = object;

    scene = unionSDF(scene, object);                                                                // Assembling scene...

    incident = light.pos - object.P;
    object2 = raymarchSDF(object2, object2.P + 2.0f*EPSILON*object2.N, normalize(incident), h);
    
    if(object2.sdf < 0.95f*length(incident))
    {
      scene.shd = 0.1f;
      //first = false;
    }
    else
    {
      scene.shd = 1.0f;
    }
  }

  // COMPUTING LIGHTNING:
  P = scene.P;
  N = scene.N;

  view = normalize(camera.pos - P);                                                                 // Computing view direction...
  incident = normalize(light.pos - P);                                                              // Computing incident light direction...
  reflected = reflect(-incident, N);                                                                // Computing reflected light direction...
  halfway = normalize(incident + view);                                                             // Computing halfway vector (Blinn-Phong)...
  ambient = light.amb;                                                                              // Computing light ambient color...
  reflection = pow(max(dot(N, halfway), 0.0f), scene.ref.w);                                        // Computing light reflection intensity

  shadow = scene.shd;
  //shadow = 1.0f;
  diffusion = clamp(dot(N, incident), 0.0f, 1.0f)*shadow;                                           // Computing light diffusion intensity... 
  
  if(i == 512 && j == 384)
  {
  printf("\n sdf = %f, inc = %f, shadow = %f\n", object2.sdf, length(incident2), shadow);
  }
  /*
  refracted = refract(incident, N, 1.00f, 1.33f);
  d = raymarch(P, -refracted, ball, 3).w;                                                           // Computing scene distance...
  P = P + d*ray;                                                                                    // Computing scene position...
  N = normal(P, ball, 3);                                                                           // Computing scene normal...
  */

  fragment.amb = ambient*light.col*scene.amb.xyz;                                                   // Computing light ambient color...
  fragment.dif = diffusion*light.col*scene.dif.xyz;                                                 // Computing light diffused color...
  fragment.ref = reflection*light.col*scene.ref.xyz;                                                // Computing light reflected color...

  //fragment_color[i + (uint)W*j] = (float4)(fragment.amb + fragment.dif, 1.0f);
  fragment_color[i + (uint)W*j] = (float4)(fragment.amb + fragment.dif + fragment.ref, 1.0f);       // Setting output color...
  //fragment_color[i + (uint)W*j] = (float4)(scene.sdf/10.0f, scene.sdf/10.0f, scene.sdf/10.0f, 1.0f);
  //fragment_color[i + (uint)W*j] = (float4)(N, 1.0f);
}
