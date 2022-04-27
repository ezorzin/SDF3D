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
  float y = 3.0f;                                                               // Sphere y-coordinate.
  float z = 0.0f;                                                               // Sphere z-coordinate.
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
  // Ray marching variables:
  int   i;                                                                      // Step index.
  float marching_distance;                                                      // Marching distance.
  vec3  marching_ray;                                                           // Marching ray.
  vec3  light_incident_ray;                                                     // Light incident ray.
  vec3  light_reflected_ray;                                                    // Light reflected ray.
  float sdf = INF;                                                              // Signed distance field.
  vec3  camera_position;                                                        // Camera position (ray marching origin).
  vec3  camera_direction;                                                       // Camera direction.
  
  // Light variables:
  vec3  light_position;                                                         // Light position.
  vec3  light_color;                                                            // Light color.
  
  float light_ambient_intensity;                                                // Light ambient intensity.
  float light_diffused_intensity;                                               // Light diffuse intensity.          
  float light_reflected_intensity;                                              // Light specular intensity.
  vec3  light_ambient_color;                                                    // Light ambient color.
  vec3  light_diffuse_color;                                                    // Light diffuse color.
  vec3  light_reflected_color;                                                  // Light reflected color.
  
  // Material variables:
  vec3  material_ambient;                                                       // Material ambient color.
  vec3  material_diffuse;                                                       // Material diffuse color.
  vec3  material_specular;                                                      // Material specular color.
  float material_shininess;                                                     // Material shininess.

  // Normal variables:
  float sdf_L = INF;                                                            // Signed distance field (left limit).
  float sdf_R = INF;                                                            // Signed distance field (right limit).
  vec3  dx = vec3(EPSILON, 0.0f, 0.0f);                                         // x-direction increment.
  vec3  dy = vec3(0.0f, EPSILON, 0.0f);                                         // y-direction increment.
  vec3  dz = vec3(0.0f, 0.0f, EPSILON);                                         // z-direction increment.
  float nx = INF;                                                               // Scene SDF's gradient x-component.
  float ny = INF;                                                               // Scene SDF's gradient y-component.
  float nz = INF;                                                               // Scene SDF's gradient z-component.
  vec3  n;                                                                      // Scene's normal vector.

  // INITIALIZING RAY MARCHING:
  camera_position = vec3(0.0f, 1.0f, 0.0f);                                     // Setting camera position (ray origin)...
  camera_direction = normalize(vec3(quad.x*AR, quad.y, 2.0f));                  // Computing position on canvas (ray intersection on quad)...
  camera_direction = normalize((inverse(V_mat)*vec4(camera_direction, 0.0f)).xyz);
  camera_position = (inverse(V_mat)*vec4(camera_position, 1.0f)).xyz;
  marching_distance = 0.0f;                                                     // Resetting marching distance...
  i = 0;                                                                        // Resetting step index... 

  // INITIALIZING LIGHT:
  light_position = vec3(5.0f, 5.0f, 0.0f);                                      // Setting light position...
  light_color = vec3(0.7f, 0.7f, 0.7f);                                         // Setting light color...
  light_ambient_intensity = 0.1f;                                               // Setting light ambient intensity...

  // INITIALIZING MATERIAL:
  material_ambient   = vec3(0.0f, 0.2f, 0.8f);                                  // Setting material ambient color...
  material_diffuse   = vec3(0.0f, 0.2f, 0.8f);                                  // Setting material diffuse color...
  material_specular  = vec3(0.5f, 0.5f, 0.5f);                                  // Setting material specular color...
  material_shininess = 32.0f;                                                   // Setting material shininess...

  // COMPUTING RAY MARCHING:
  while (
          (marching_distance < MAX_DISTANCE) &&                                 // Checking matching distance...
          (sdf > EPSILON) &&                                                    // Checking signed distance field...
          (i < MAX_STEPS)                                                       // Checking step index...
        )
  {
    marching_ray = camera_position + marching_distance*camera_direction;        // Computing marching ray...      
    sdf = INF;                                                                  // Setting SDF to infinity...
    sdf = min(sdf, SDF_plane(marching_ray));                                    // Computing signed distance field...
    sdf = min(sdf, SDF_sphere(marching_ray));                                   // Computing signed distance field...
    marching_distance += sdf;                                                   // Updating marching distance...
    i++;                                                                        // Updating step index...
  }

  marching_ray = camera_position + marching_distance*camera_direction;          // Computing final marching ray...

  // COMPUTING NORMAL VECTOR:
  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(marching_ray - dx));                             // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(marching_ray - dx));                            // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(marching_ray + dx));                             // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(marching_ray + dx));                            // Computing signed distance field (right limit)...
  
  nx = sdf_R - sdf_L;                                                           // Computing gradient (x-component)...      
  
  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(marching_ray - dy));                             // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(marching_ray - dy));                            // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(marching_ray + dy));                             // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(marching_ray + dy));                            // Computing signed distance field (right limit)...
  
  ny = sdf_R - sdf_L;                                                           // Computing gradient (y-component)... 

  sdf_L = INF;                                                                  // Setting SDF to infinity (left limit)...
  sdf_L = min(sdf_L, SDF_plane(marching_ray - dz));                             // Computing signed distance field (left limit)...
  sdf_L = min(sdf_L, SDF_sphere(marching_ray - dz));                            // Computing signed distance field (left limit)...
  
  sdf_R = INF;                                                                  // Setting SDF to infinity (right limit)...
  sdf_R = min(sdf_R, SDF_plane(marching_ray + dz));                             // Computing signed distance field (right limit)...
  sdf_R = min(sdf_R, SDF_sphere(marching_ray + dz));                            // Computing signed distance field (right limit)...
  
  nz = sdf_R - sdf_L;                                                           // Computing gradient (z-component)... 
  
  n = normalize(vec3(nx, ny, nz));                                              // Computing normal vector...
  
  // COMPUTING LIGHTNING:
  light_incident_ray = normalize(light_position - marching_ray);                // Computing light direction...
  light_reflected_ray = reflect(-light_incident_ray, n);                      // Computing light reflection direction...
  
  light_ambient_color = light_ambient_intensity*light_color;
  light_diffused_intensity = clamp(dot(n, light_incident_ray), 0.0f, 1.0f);         // Computing light diffusion intensity...
  light_diffuse_color = light_color*(light_diffused_intensity*material_diffuse); // Computing diffuse light color...
  
  //light_reflected_intensity = pow(max(dot(camera_position - marching_ray, light_reflection_direction), 0.0), material_shininess);
  //light_reflected_intensity = pow(max(dot(-camera_direction, light_reflection_direction), 0.0), material_shininess);
  light_reflected_intensity = pow(max(dot(normalize(camera_position - marching_ray), light_reflected_ray), 0.0), material_shininess);
  light_reflected_color = light_color*(light_reflected_intensity*material_specular);    // Computing specular light color...

  fragment_color = vec4(light_ambient_color + light_diffuse_color + light_reflected_color, 1.0f);  // Setting output color...
}