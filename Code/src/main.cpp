/// @file

#define INTEROP       true                                                                          // "true" = use OpenGL-OpenCL interoperability.
#define SX            800                                                                           // Window x-size [px].
#define SY            600                                                                           // Window y-size [px].
#define NM            "Neutrino - 3D signed distance field"                                         // Window name.
#define OX            0.0f                                                                          // x-axis orbit initial rotation.
#define OY            0.0f                                                                          // y-axis orbit initial rotation.
#define PX            0.0f                                                                          // x-axis pan initial translation.
#define PY            0.0f                                                                          // y-axis pan initial translation.
#define PZ            0.0f                                                                          // z-axis pan initial translation.

#ifdef __linux__
  #define SHADER_HOME "../../Code/shader/"                                                          // Linux OpenGL shaders directory.
  #define KERNEL_HOME "../../Code/kernel/"                                                          // Linux OpenCL kernels directory.
  #define GMSH_HOME   "../../Code/mesh/"                                                            // Linux GMSH mesh directory.
#endif

#ifdef WIN32
  #define SHADER_HOME "..\\..\\Code\\shader\\"                                                      // Windows OpenGL shaders directory.
  #define KERNEL_HOME "..\\..\\Code\\kernel\\"                                                      // Windows OpenCL kernels directory.
  #define GMSH_HOME   "..\\..\\Code\\mesh\\"                                                        // Windows GMSH mesh directory.
#endif

#define SHADER_VERT   "voxel_vertex.vert"                                                           // OpenGL vertex shader.
#define SHADER_GEOM   "voxel_geometry.geom"                                                         // OpenGL geometry shader.
#define SHADER_FRAG   "voxel_fragment.frag"                                                         // OpenGL fragment shader.
#define KERNEL_1      "thekernel_1.cl"                                                              // OpenCL kernel source.
#define UTILITIES     "utilities.cl"                                                                // OpenCL utilities source.

// INCLUDES:
#include "nu.hpp"                                                                                   // Neutrino's header file.

int main ()
{
  // MOUSE PARAMETERS:
  float               ms_orbit_rate  = 1.0f;                                                        // Orbit rotation rate [rev/s].
  float               ms_pan_rate    = 5.0f;                                                        // Pan translation rate [m/s].
  float               ms_decaytime   = 1.25f;                                                       // Pan LP filter decay time [s].

  // GAMEPAD PARAMETERS:
  float               gmp_orbit_rate = 1.0f;                                                        // Orbit angular rate coefficient [rev/s].
  float               gmp_pan_rate   = 1.0f;                                                        // Pan translation rate [m/s].
  float               gmp_decaytime  = 1.25f;                                                       // Low pass filter decay time [s].
  float               gmp_deadzone   = 0.30f;                                                       // Gamepad joystick deadzone [0...1].

  // OPENGL:
  nu::opengl*         gl             = new nu::opengl (NM, SX, SY, OX, OY, PX, PY, PZ);             // OpenGL context.
  nu::shader*         sh             = new nu::shader ();                                           // OpenGL shader program.
  nu::projection_mode pmode          = nu::MONOCULAR;                                               // OpenGL projection mode.
  nu::view_mode       vmode          = nu::INVERSE;                                                 // OpenGL VIEW mode.

  // OPENCL:
  nu::opencl*         cl             = new nu::opencl (nu::GPU);                                    // OpenCL context.
  nu::kernel*         K1             = new nu::kernel ();                                           // OpenCL kernel array.
  nu::float4*         color          = new nu::float4 (0);                                          // Color [].
  nu::float4*         V              = new nu::float4 (1);                                          // View matrix.
  nu::float4*         canvas_param   = new nu::float4 (2);                                          // Canvas parameters.
  nu::float4*         camera_param   = new nu::float4 (3);                                          // Camera parameters.
  nu::float4*         light_param    = new nu::float4 (4);                                          // Light parameters.
  nu::float4*         material_param = new nu::float4 (5);                                          // Material parameters.

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// DATA INITIALIZATION //////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  int                 i;

  for(i = 0; i < (800*600); i++)
  {
    color->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                               // Setting node color...
  }

  V->data.push_back ({1.0f, 0.0f, 0.0f, 0.0f});
  V->data.push_back ({0.0f, 1.0f, 0.0f, 0.0f});
  V->data.push_back ({0.0f, 0.0f, 1.0f, 0.0f});
  V->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});

  canvas_param->data.push_back ({800.0f, 600.0f, 800.0f/600.0f, 60.0f});
  camera_param->data.push_back ({0.0f, 0.2f, 2.0f, 1.0f});
  light_param->data.push_back ({5.0f, 5.0f, 0.0f, 1.0f});
  light_param->data.push_back ({0.7f, 0.7f, 0.7f, 0.1f});
  material_param->data.push_back ({0.0f, 0.2f, 0.8f, 1.0f});
  material_param->data.push_back ({0.0f, 0.2f, 0.8f, 1.0f});
  material_param->data.push_back ({0.5f, 0.5f, 0.5f, 1.0f});
  material_param->data.push_back ({12.0f, 12.0f, 12.0f, 12.0f});

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// OPENCL KERNELS INITIALIZATION /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  K1->addsource (std::string (KERNEL_HOME) + std::string (UTILITIES));                              // Setting kernel source file...
  K1->addsource (std::string (KERNEL_HOME) + std::string (KERNEL_1));                               // Setting kernel source file...
  K1->build (SX, SY, 0);                                                                            // Building kernel program...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// OPENGL SHADERS INITIALIZATION /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  sh->addsource (std::string (SHADER_HOME) + std::string (SHADER_VERT), nu::VERTEX);                // Setting shader source file...
  sh->addsource (std::string (SHADER_HOME) + std::string (SHADER_GEOM), nu::GEOMETRY);              // Setting shader source file...
  sh->addsource (std::string (SHADER_HOME) + std::string (SHADER_FRAG), nu::FRAGMENT);              // Setting shader source file...
  sh->build (1);                                                                                    // Building shader program...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// SETTING OPENCL KERNEL ARGUMENTS /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  cl->write ();                                                                                     // Writing OpenCL data...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// APPLICATION LOOP ////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  while(!gl->closed ())                                                                             // Opening window...
  {
    cl->get_tic ();                                                                                 // Getting "tic" [us]...

    V->data[0] = {gl->V_mat[0], gl->V_mat[4], gl->V_mat[8], gl->V_mat[12]};
    V->data[1] = {gl->V_mat[1], gl->V_mat[5], gl->V_mat[9], gl->V_mat[13]};
    V->data[2] = {gl->V_mat[2], gl->V_mat[6], gl->V_mat[10], gl->V_mat[14]};
    V->data[3] = {gl->V_mat[3], gl->V_mat[7], gl->V_mat[11], gl->V_mat[15]};

    //std::cout << "T0 = " << gl->V_mat[12] << "T1 = " << gl->V_mat[13] << "T2 = " << gl->V_mat[14] <<
    //  "T3 = " << gl->V_mat[15] << std::endl;

    cl->write (1);
    cl->acquire ();                                                                                 // Acquiring OpenCL kernel...
    cl->execute (K1, nu::WAIT);                                                                     // Executing OpenCL kernel...
    cl->release ();                                                                                 // Releasing OpenCL kernel...

    gl->begin ();                                                                                   // Beginning gl...
    gl->poll_events ();                                                                             // Polling gl events...
    gl->mouse_navigation (ms_orbit_rate, ms_pan_rate, ms_decaytime);                                // Polling mouse...
    gl->gamepad_navigation (gmp_orbit_rate, gmp_pan_rate, gmp_decaytime, gmp_deadzone);             // Polling gamepad...
    gl->plot (sh, pmode, vmode);                                                                    // Plotting shared arguments...
    gl->end ();                                                                                     // Ending gl...

    cl->get_toc ();                                                                                 // Getting "toc" [us]...
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////// CLEANUP ////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  delete cl;                                                                                        // Deleting OpenCL context...
  delete gl;                                                                                        // Deleting OpenGL context...
  delete sh;                                                                                        // Deleting shader...
  delete K1;
  delete color;                                                                                     // Deleting color data...
  delete color;
  delete V;
  delete camera_param;
  delete canvas_param;
  delete light_param;
  delete material_param;

  return 0;
}
