/// @file

#define INTEROP       true                                                                          // "true" = use OpenGL-OpenCL interoperability.
#define SX            1024                                                                          // Window x-size [px].
#define SY            768                                                                           // Window y-size [px].
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
  nu::float4*         fragment_color = new nu::float4 (0);                                          // Color [r, g, b, a].
  nu::float4*         center         = new nu::float4 (1);                                          // Center point.
  nu::float16*        view_matrix    = new nu::float16 (2);                                         // View matrix [4x4].
  nu::float4*         canvas         = new nu::float4 (3);                                          // Canvas parameters [W, H, AR, FOV].
  nu::float4*         light_position = new nu::float4 (4);                                          // Light position [x, y, z, k].
  nu::float4*         light_color    = new nu::float4 (5);                                          // Light color [r, g, b, ambient].
  nu::int1*           object_number  = new nu::int1 (6);                                            // Number of objects.
  nu::int1*           object_type    = new nu::int1 (7);                                            // Object type.
  nu::float16*        T              = new nu::float16 (8);                                         // Transformation matrix.
  nu::float16*        M              = new nu::float16 (9);                                         // Material matrix.

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// DATA INITIALIZATION //////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  int                 i;

  for(i = 0; i < (SX*SY); i++)
  {
    fragment_color->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                      // Initializing fragment color...
  }

  center->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});

  view_matrix->data.push_back (
                               {1.0f, 0.0f, 0.0f, 0.0f,                                             // Initializing view matrix...
                                0.0f, 1.0f, 0.0f, 0.0f,                                             // Initializing view matrix...
                                0.0f, 0.0f, 1.0f, 0.0f,                                             // Initializing view matrix...
                                0.0f, 0.0f, 0.0f, 1.0f}                                             // Initializing view matrix...
                              );

  canvas->data.push_back ({800.0f, 800.0f, 800.0f/600.0f, 60.0f});                                  // Initializing canvas parameters [W, H, AR, FOV]...
  light_position->data.push_back ({5.0f, 5.0f, 0.0f, 10.0f});                                       // Initializing light position [x, y, z, k]...
  light_color->data.push_back ({1.0f, 1.0f, 1.0f, 0.1f});                                           // Initializing light color [r, g, b, ambient]...
  object_number->data.push_back (3);

  object_type->data.push_back (2);
  object_type->data.push_back (2);
  object_type->data.push_back (2);

  for(i = 0; i < 3; i++)
  {
    T->data.push_back (
                       {1.0f, 0.0f, 0.0f, float(i),
                        0.0f, 1.0f, 0.0f, 0.0f,
                        0.0f, 0.0f, 1.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 1.0f}
                      );
    M->data.push_back (
                       {0.2f, 0.0f, 0.0f, 0.0f,
                        0.0f, 0.8f, 0.2f, 1.0f,
                        0.0f, 0.8f, 0.2f, 1.0f,
                        0.5f, 0.5f, 0.5f, 12.0f}
                      );
  }

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

    view_matrix->data[0] = {gl->V_mat[0], gl->V_mat[4], gl->V_mat[8],  gl->V_mat[12],
                            gl->V_mat[1], gl->V_mat[5], gl->V_mat[9],  gl->V_mat[13],
                            gl->V_mat[2], gl->V_mat[6], gl->V_mat[10], gl->V_mat[14],
                            gl->V_mat[3], gl->V_mat[7], gl->V_mat[11], gl->V_mat[15]};
    canvas->data[0]      = {float(SX), float(SY), gl->aspect_ratio, 60.0f};

    cl->write (2);
    cl->write (3);

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
  delete fragment_color;
  delete center;
  delete view_matrix;
  delete canvas;
  delete light_position;
  delete light_color;                                                                               // Deleting color data...
  delete object_number;
  delete object_type;
  delete T;
  delete M;

  return 0;
}
