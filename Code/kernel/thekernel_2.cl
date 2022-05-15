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

/////////////////////////////////////////////////////////////////////////////////////////////////////
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
  uint            k = get_global_id(2);                                                             // Global index k [#].
  float16         V;                                                                                // View matrix.
  float           W = canvas[0].x;                                                                  // Window width [px].
  float           H = canvas[0].y;                                                                  // Window height [px].
  float           AR = canvas[0].z;                                                                 // Window aspect ratio [].
  float           FOV = canvas[0].w;                                                                // FOV [degrees].
  float           x = 2.0f*(i/W) - 1.0f;                                                            // Canvas x-coordinate [-1.0f...+1.0f].
  float           y = 2.0f*(j/H) - 1.0f;                                                            // Canvas y-coordinate [-1.0f...+1.0f].
  int             n = object_number[0];

  int             m;
  float3          P;                                                                                // Scene position.
  float3          N;                                                                                // Scene normal direction.

  struct Camera   camera;                                                                           // Camera.
  struct Light    light;                                                                            // Light.
  struct Object   scene;                                                                            // Scene.
  struct Fragment fragment;                                                                         // Fragment.
  struct Object   object;


  // INITIALIZING RAY MARCHING:
  V = view_matrix[0];                                                                               // Getting view matrix...

  camera.pos = (float3)(0.0f, 0.2f, 2.0f);                                                          // Setting initial camera position...
  camera.pos = mul(V, (float4)(camera.pos, 1.0f)).xyz;                                              // Applying arcball to camera position...
  camera.fov = FOV;                                                                                 // Setting camera FOV...

  light.pos = light_position[0].xyz;                                                                // Setting light position...
  light.col = light_color[0].xyz;                                                                   // Setting light color [r, g, b]...
  light.amb = light_color[0].w;                                                                     // Setting light ambient intensity [ambient]...
  light.k   = light_position[0].w;                                                                  // Setting light sharpness [k]...

  ray = (float3)(x*AR, y, -1.0f/tan(camera.fov*PI/360.0f));                                         // Computing ray intersection on canvas...
  ray = mul(V, (float4)(ray, 1.0f)).xyz;                                                            // Applying arcball to ray direction...
  direction = normalize(ray);                                                                             // Normalizing ray direction...

  P = camera.pos + distance*direction;

  for (m = 0; m < object_number; m++)
  {
    object.type = object_type[k];                                                                   // Getting object type...
    object.T = T[k];                                                                                // Getting object transformation matrix...
    object.par = M[k].s0123;                                                                        // Getting object parameters...
    object.amb = M[k].s4567;                                                                        // Getting object ambient color...
    object.dif = M[k].s89AB;                                                                        // Getting object diffusion color...
    object.ref = M[k].sCDEF;                                                                        // Getting object reflection color...
    
    object.sdf = objectSDF(object.type, object.par, object.T, P);                                                             // Computing object sdf...
    scene = unionSDF(scene, object);                                                                // Assembling scene...

    object_dSDF.x = objectSDF(object.type, object.par, object.T, P + TET_0*h);                                                     // Computing sdf (X left differential)...
    object_dSDF.y = objectSDF(object.type, object.par, object.T, P + TET_1*h);                                                     // Computing sdf (X right differential)...
    object_dSDF.z = objectSDF(object.type, object.par, object.T, P + TET_2*h);                                                     // Computing sdf (Y left differential)...
    object_dSDF.w = objectSDF(object.type, object.par, object.T, P + TET_3*h);                                                     // Computing sdf (Y right differential)...
    
    scene_dSDF.x = unionSDFnorm(scene_dSDF.x, object_dSDF.x);                                                      
    scene_dSDF.y = unionSDFnorm(scene_dSDF.y, object_dSDF.y);
    scene_dSDF.z = unionSDFnorm(scene_dSDF.z, object_dSDF.z);
    scene_dSDF.w = unionSDFnorm(scene_dSDF.w, object_dSDF.w);
  }

}


