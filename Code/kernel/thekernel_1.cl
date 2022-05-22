/// @file
// CONSTANTS:
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

#define W canvas[0].x                                                                               // Window width [px].
#define H canvas[0].y                                                                               // Window height [px].
#define AR canvas[0].z                                                                              // Window aspect ratio [].
#define FOV canvas[0].w                                                                             // FOV [degrees].
#define X 2.0f*(get_global_id(0)/W) - 1.0f                                                          // Canvas x-coordinate [-1.0f...+1.0f].
#define Y 2.0f*(get_global_id(1)/H) - 1.0f                                                          // Canvas y-coordinate [-1.0f...+1.0f].

#define V_0123 V[0].s0123
#define V_4567 V[0].s4567
#define V_89AB V[0].s89AB
#define V_CDEF V[0].sCDEF

#define CAMERA_I (float4)(0.0f, 0.2f, 2.0f, 1.0f)                                                   // Initial camera position.
#define CAMERA_POS_X dot(V_0123, CAMERA_I)
#define CAMERA_POS_Y dot(V_4567, CAMERA_I)
#define CAMERA_POS_Z dot(V_89AB, CAMERA_I)
#define CAMERA_POS (float3)(CAMERA_POS_X, CAMERA_POS_Y, CAMERA_POS_Z)                               // Arcball applied to camera position.

#define LIGHT_POS light_position[0].xyz                                                             // Light position.
#define LIGHT_K light_position[0].w                                                                 // Light sharpness [k]...
#define LIGHT_COL light_color[0].xyz                                                                // Light color [r, g, b]...
#define LIGHT_AMB light_color[0].w                                                                  // Light ambient intensity [ambient]...

#define PAR_0 M[k].s0
#define PAR_1 M[k].s1
#define PAR_2 M[k].s2
#define PAR_3 M[k].s3
#define AMB_R M[k].s4
#define AMB_G M[k].s5
#define AMB_B M[k].s6
#define ALPHA M[k].s7
#define DIF_R M[k].s8
#define DIF_G M[k].s9
#define DIF_B M[k].sA
#define N_RATIO M[k].sB
#define REF_R M[k].sC
#define REF_G M[k].sD
#define REF_B M[k].sE
#define K M.sF

#define AMB M[k].s4567
#define DIF M[k].s89AB
#define REF M[k].sCDEF

#define T_0123 T[k].s0123
#define T_4567 T[k].s4567
#define T_89AB T[k].s89AB
#define T_CDEF T[k].sCDEF

#define RAY_TMP float3_tmp       
#define P_TMP float3_tmp             

#define TYPE object_type[k]
#define P object_position[k]
#define N object_normal[k]
#define SDF object_sdf[k]

#define SCENE 0
#define PLANE 1
#define SPHERE 2

#define SDF_SCENE INF
#define SDF_PLANE P.y
#define SDF_SPHERE (length(P) - PAR_0)

#define N_SCENE ZERO3
#define N_PLANE ((P + TET_A*EPSILON).y*TET_A + (P + TET_B*EPSILON).y*TET_B + (P + TET_C*EPSILON).y*TET_C + (P + TET_D*EPSILON).y*TET_D)
#define N_SPHERE ((length(P + TET_A*EPSILON) - PAR_0)*TET_A + (length(P + TET_B*EPSILON) - PAR_0)*TET_B + (length(P + TET_C*EPSILON) - PAR_0)*TET_D + (length(P + TET_D*EPSILON) - PAR_0)*TET_D)



//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////// KERNEL ///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
__kernel void thekernel(__global float4*    fragment_color,                                         // Fragment color.
                        __global float4*    center,
                        __global float16*   V,                                                      // View matrix [4x4].
                        __global float4*    canvas,                                                 // Canvas [W, H, AR, FOV].
                        __global float4*    light_position,                                         // Light position [x, y, z, k].
                        __global float4*    light_color,                                            // Light color [r, g, b, ambient].
                        __global int*       object_type,                                            // Object type.
                        __global float3*    object_position,
                        __global float3*    object_normal,
                        __global float*     object_sdf,                                             // Object sdf.
                        __global float16*   T,                                                      // Object transformation matrix [4x4].                     
                        __global float16*   M                                                       // Object material matrix [4x4].
                        )
{ 
  uint            k = get_global_id(2);                                                             // Global index j [#].
 
  uint            n;                                                                                // Marching index.
  float3          float3_tmp;
  float3          ray;                                                                              // Marching ray.
  float           distance;                                                                         // Marching distance.
  
  float3 p;
  uint q;

  // INITIALIZING RAY MARCHING:
  ray = (float3)(X*AR, Y, -1.0f/tan(FOV*M_PI_F/360.0f));                                            // Computing ray intersection on canvas...
  RAY_TMP.x = dot(V_0123, (float4)(ray, 1.0f));                                                     // Applying arcball to ray direction...
  RAY_TMP.y = dot(V_4567, (float4)(ray, 1.0f));                                                     // Applying arcball to ray direction...
  RAY_TMP.z = dot(V_89AB, (float4)(ray, 1.0f));                                                     // Applying arcball to ray direction...                                                          
  ray = normalize(RAY_TMP);                                                                         // Normalizing ray direction...

  // COMPUTING RAY MARCHING:
  distance = 0.0f;                                                                                  // Marching distance.

  for (q = 0; q < MAX_STEPS; q++)
  {
  for (n = 0; n < MAX_STEPS; n++)
  {
    p = CAMERA_POS + distance*ray;                                                                  // Computing marching ray...
    P_TMP.x = dot(T_0123, (float4)(p, 1.0f));                                                       // Applying transformation matrix...
    P_TMP.y = dot(T_4567, (float4)(p, 1.0f));                                                       // Applying transformation matrix...
    P_TMP.z = dot(T_89AB, (float4)(p, 1.0f));                                                       // Applying transformation matrix...                                                          
    p = P_TMP;                                                                                      // Applying transformation matrix... 
  }

  object_position[k] = p;
  }
  /*
    switch (TYPE)
    {
      case SCENE:
        SDF = SDF_SCENE;
        break;

      case PLANE:
        SDF = SDF_PLANE;
        break;

      case SPHERE:
        SDF = SDF_SPHERE;                                                         
        break;
    }
 
    distance += SDF;                                                                      // Updating marching distance...
    
    if(distance > MAX_DISTANCE || SDF < EPSILON) break;                                   // Checking numeric precision constraints...
  }*/
  /*
  if(distance > MAX_DISTANCE)
  {
    AMB = ZERO4;
    DIF = ZERO4;
    REF = ZERO4;
  }

  SDF = distance;                                                                            // Computing object distance...
  P = CAMERA_POS + distance*ray;                                                                  // Computing object position...                                          
  
  switch (TYPE)
  {
    case SCENE:
      N = N_SCENE;

    case PLANE:                                        
      N = N_PLANE;
      break;

    case SPHERE:
      N = N_SPHERE;
      break;
  }
  */
}
