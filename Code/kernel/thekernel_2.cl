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

  struct Object unionSDF(struct Object a, struct Object b)
{
  return a.sdf < b.sdf? a : b;
}

struct Object smoothunionSDF(struct Object A, struct Object B, float smoothness)
{
  float interpolation = clamp(0.5f + 0.5f*(B.sdf - A.sdf)/smoothness, 0.0f, 1.0f);
  struct Object C;

  C = A.sdf < B.sdf? A : B;

  
  C.amb = mix(B.amb, A.amb, interpolation);
  C.dif = mix(B.dif, A.dif, interpolation);
  C.ref = mix(B.ref, A.ref, interpolation);
  //C.P = mix(B.P, A.P, interpolation);
  //C.N = mix(B.N, A.N, interpolation);
  
  //C.sdf = mix(B.sdf, A.sdf, interpolation);


  return C;
}

struct Object shadowSDF(struct Object a, struct Object b)
{
  return a.shd < b.shd? a : b;
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
/////////////////////////////////////////// SHADOW MARCHING //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
struct Object shadowmarchSDF(struct Object object, float3 origin, float3 direction, float h, struct Light light)
{
  float             distance = 0.0f;                                                                // Marching distance.
  float             sdf_prev = INF;                                                                 // Previous step signed distance field.
  int               i;                                                                              // Step index.
  float3            P;                                                                              // Marching ray.
  float             intersection;                                                                   // Previous to current SDF on-ray intersection point.
  float             distance_est;                                                                   // Estimated closest distance.

  object.shd = 1.0f;                                                                                // Shadow intensity.

  for (i = 0; i < MAX_STEPS; i++)
  {
    P = origin + distance*direction;                                                                // Computing marching ray...   
    object.sdf = objectSDF(object, P);                                                              // Computing scene sdf...                                          
    intersection = (i == 0) ? 0.0f : pown(object.sdf, 2)/(2.0f*sdf_prev);                           // Computing on-ray intersection point...
    distance_est = sqrt(pown(object.sdf, 2) - pown(intersection, 2));                               // Computing estimated closest distance...
    object.shd = min(object.shd, light.k*distance_est/max(0.0f, distance - intersection));          // Computing shadow intensity...
    sdf_prev = object.sdf;                                                                          // Backing up signed distance field...
    distance += object.sdf;                                                                         // Updating marching distance...
    
    if(distance > MAX_DISTANCE || object.shd < EPSILON) break;                                      // Checking numeric precision constraints...
  }

  object.shd = clamp(object.shd, 0.0f, 1.0f);                                                       // Clamping shadow intensity...
  object.sdf = distance;

  return object;                                                                                    // Returning final shadow intensity...
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

scene = unionSDF(scene, object);                                                                // Assembling scene... EZOR: to be done in a second kernel.
  

  // COMPUTING LIGHTNING:
  P = scene.P;
  N = scene.N;

  view = normalize(CAMERA_POS - P);                                                                 // Computing view direction...
  incident = normalize(light.pos - P);                                                              // Computing incident light direction...
  reflected = reflect(-incident, N);                                                                // Computing reflected light direction...
  halfway = normalize(incident + view);                                                             // Computing halfway vector (Blinn-Phong)...
  ambient = light.amb;                                                                              // Computing light ambient color...
  reflection = pow(max(dot(N, halfway), 0.0f), scene.ref.w);                                        // Computing light reflection intensity

  scene2.shd = 1.0f;

  /*
  // COMPUTING SHADOW MARCHING:
  for(k = 0; k < n; k++)
  {
    object.type = object_type[k];                                                                   // Getting object type...
    object.T = T[k];                                                                                // Getting object transformation matrix...
    object.par = M[k].s0123;                                                                        // Getting object parameters...
    object = shadowmarchSDF(object, P + 2.0f*EPSILON*N, incident, h, light);
    scene2 = shadowSDF(scene2, object);
  }

  shadow = scene2.shd;
  */
  shadow = 1.0f;

  diffusion = clamp(dot(N, incident), 0.0f, 1.0f)*shadow;                                           // Computing light diffusion intensity... 
  
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