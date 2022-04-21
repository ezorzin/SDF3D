/// @file

#version 460 core

uniform mat4 V_mat;                                                             // View matrix.
uniform mat4 P_mat;                                                             // Projection matrix.
uniform float AR;                                                               // Framebuffer aspect ratio.

in  vec4 color;                                                                 // Voxel color.
in  vec2 quad;                                                                  // Voxel quad.

out vec4 fragment_color;                                                        // Fragment color.
 
// CONSTANTS:
#define PI 3.1415925359f                                                        // PI.
#define TWO_PI 6.2831852f                                                       // 2*PI.
#define MAX_STEPS 100                                                           // Maximum ray marching steps.
#define MAX_DISTANCE 100.0f                                                     // Maximum marching distance.
#define EPSILON 0.01f                                                           // Minimum surface distance (ray marching epsilon).
#define INF 1.0f/0.0f                                                           // Infinity.
 
//////////////////////////////////////////////////////////////////////////////////
/////////////// OBJECT'S SIGNED DISTANCE FIELDS IMPLICIT FUNCTIONS ///////////////
//////////////////////////////////////////////////////////////////////////////////
float SDF_sphere(vec3 p)
{
  float x = 0.0f;                                                               // Sphere x-coordinate.
  float y = 1.0f;                                                               // Sphere y-coordinate.
  float z = 6.0f;                                                               // Sphere z-coordinate.
  float r = 1.0f;                                                               // Sphere radius.
  vec3  s = vec3(x, y, z);                                                      // Sphere position.
  float SDF = length(p - s) - r;                                                // Signed distance field.
 
  return SDF;                                                                   // Returning signed distance field...
}

float SDF_plane(vec3 p)
{
  float SDF = p.y;                                                              // Signed distance field.
 
  return SDF;                                                                   // Returning signed distance field...
}

//////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// SCENE RENDERING ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void main()
{
  int   i;                                                                      // Step index.
  float d_marching;                                                             // Marching distance.
  vec3  ray_marching;                                                           // Marching ray.
  float sdf = INF;                                                              // Signed distance field.
  float sdf_L = INF;                                                            // Signed distance field (left limit).
  float sdf_R = INF;                                                            // Signed distance field (right limit).
  vec3  ray_origin;                                                             // Ray origin (vantage point).
  vec3  ray_direction;                                                          // Ray direction.
  vec3  light_position;                                                         // Light position.
  vec3  light_direction;                                                        // Light direction.
  float light_intensity;                                                        // Light intensity.
  float light_ambient;                                                          // Ambient light intenstity.
  float light_diffuse;                                                          // Diffuse light intensity.
  float light_specular;                                                         // Specular light intensity.
  vec3  material_ambient;
  vec3  material_diffuse;
  vec3  material_specular;
  float material_shininess;
  vec3  dx = vec3(EPSILON, 0.0f, 0.0f);                                         // x-direction increment.
  vec3  dy = vec3(0.0f, EPSILON, 0.0f);                                         // y-direction increment.
  vec3  dz = vec3(0.0f, 0.0f, EPSILON);                                         // z-direction increment.
  float nx = INF;                                                               // Scene SDF's gradient x-component.
  float ny = INF;                                                               // Scene SDF's gradient y-component.
  float nz = INF;                                                               // Scene SDF's gradient z-component.
  vec3  n;                                                                      // Scene's normal vector.

  // INITIALIZING SCENE:
  light_position = vec3(5.0f, 5.0f, 0.0f);                                      // Setting light position...
  ray_origin = vec3(0.0f, 1.0f, 0.0f);                                          // Setting ray origin (vantage point)...
  ray_direction = normalize(vec3(quad.x*AR, quad.y, 1.0f));                     // Computing ray direction...
  d_marching = 0.0f;                                                            // Resetting marching distance...
  i = 0;                                                                        // Resetting step index... 

  material_ambient   = vec3(1.0f, 0.5f, 0.31f);
  material_diffuse   = vec3(1.0f, 0.5f, 0.31f);
  material_specular  = vec3(0.5f, 0.5f, 0.5f);
  material_shininess = 32.0f;

  // RAY MARCHING:
  while (
          (d_marching < MAX_DISTANCE) && 
          (sdf > EPSILON) &&
          (i < MAX_STEPS)
        )
  {
    ray_marching = ray_origin + d_marching*ray_direction;                       // Computing marching ray...      
    sdf = INF;                                                                  // Setting SDF to infinity...
    sdf = min(sdf, SDF_plane(ray_marching));                                    // Computing signed distance field...
    sdf = min(sdf, SDF_sphere(ray_marching));                                   // Computing signed distance field...
    d_marching += sdf;                                                          // Updating marching distance...
    i++;
  }

  ray_marching = ray_origin + d_marching*ray_direction;                         // Computing final marching ray...

  // COMPUTING NORMAL VECTOR:
  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(ray_marching - dx));                             // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(ray_marching - dx));                            // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(ray_marching + dx));                             // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(ray_marching + dx));                            // Computing signed distance field (right limit)...
  
  nx = sdf_R - sdf_L;                                                           // Computing gradient (x-component)...      
  
  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(ray_marching - dy));                             // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(ray_marching - dy));                            // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(ray_marching + dy));                             // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(ray_marching + dy));                            // Computing signed distance field (right limit)...
  
  ny = sdf_R - sdf_L;                                                           // Computing gradient (y-component)... 

  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(ray_marching - dz));                             // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(ray_marching - dz));                            // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(ray_marching + dz));                             // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(ray_marching + dz));                            // Computing signed distance field (right limit)...
  
  nz = sdf_R - sdf_L;                                                           // Computing gradient (z-component)... 
  
  n = normalize(vec3(nx, ny, nz));                                              // Computing normal vector...
  
  // COMPUTING LIGHTNING:
  light_ambient = 0.01f;
  light_direction = normalize(light_position - ray_marching);                   // Computing light direction...
  light_diffuse = clamp(dot(n, light_direction), 0.0f, 1.0f);                         // Computing diffuse light intensity...
  
  vec3 reflectDir = reflect(-light_direction, n);
  float shininess = 50;
  float specularStrength = 0.2;
  float spec = pow(max(dot(ray_origin - ray_direction, reflectDir), 0.0), shininess);
  light_specular = specularStrength * spec; 

  light_intensity = light_ambient + light_diffuse + light_specular;

  fragment_color = vec4(light_intensity, light_intensity, light_intensity, 1.0f);                                      // Setting output color...
}