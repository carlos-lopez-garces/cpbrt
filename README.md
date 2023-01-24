**CPBRT** is my physically-based, offline toy renderer. It is the result of studying Physically Based Rendering: From Theory to Implementation, by Pharr, Jakob, and Humphreys ([online edition](https://pbrt.org/)), writing down the code fragments provided by the authors, and filling in the gaps.

A more detailed description is found in [my blog](https://carlos-lopez-garces.github.io/projects/cpbrt/).

## Features

<table>
  <tr>
    <td> <b>Light transport algorithms:</b> Kajiya path tracing (unidirectional, unbiased Monte Carlo estimation of the light transport equation). Direct-lighting (no indirect illumination) and path (full global illumination) integrators. </td>
    <td> <b>Reflectance models and BRDFs:</b> Lambert diffuse model, Oren-Nayar diffuse model for rough surfaces, Fresnel perfectly specular model, and Fresnel glossy specular model (with Torrance-Sparrow microfacets with Beckmann-Spizzichino or Trowbridge-Reitz distributions). </td>
  </tr>
  <tr>
    <td> <b>Textures:</b> Floating-point and spectrum constant-value textures. Procedural checkerboard texture, antialiased with a box filter. Mipmapping. </td>
    <td> <b>Materials:</b> Matte with either a perfect diffuse Lambertian BRDF or an Oren-Nayar BRDF for various degrees of roughness; plastic with diffuse and glossy specular BRDFs; mirror with a perfectly-specular BRDF; gold; glass with perfectly-specular BRDF and BTDF; diffuse substrate and glossy coat with an Ashikhmin-Shirley BRDF. </td>
  </tr>
  <tr>
    <td> <b>Shapes:</b> Triangle meshes, single triangles, and spherical implicit surfaces. </td>
    <td> <b>Accelerators:</b> BVH with 5 different primitive (or object) subdivision methods: linear BVH, hierarchical linear BVH, midpoint partitioning, equal counts partitioning, and surface area heuristic (SAH).</td>
  </tr>
    <td> <b>Samplers:</b> Uniform or jittered stratified pixel sampling for 1D samples and Latin Hypercube sampling for 2D samples. Samplers rely on a Permuted Congruential Generator (PCG) pseudo-random number generator. </td>
    <td> <b>Filters:</b> Box, triangle, Gaussian, Mitchell-Netravali, and Lanczos windowed-sinc filters. </td>
  <tr>
    <td> <b>Lights:</b> Point, distant, and diffuse area light sources. An area light can take the form of any of the supported *shapes*. Infinite area light source backed by environment map. </td>
    <td> <b>Cameras:</b> Thin lens perspective and orthographic projective cameras with configurable aperture and focal distance (for depth of field) and film aspect ratio. The perspective camera also has a configurable field of view. </td>
  </tr>
  <tr>
    <td> <b>Participating media:</b> Homogeneous-density and grid-based variable-density media. </td>
    <td> </td>
  </tr>
</table>

## Select images

<table>
  <tr>
    <td><p align="center">
      <i>Diffuse substrate and glossy coat with an Ashikhmin-Shirley BRDF. Volumetric path integrator.</i>
      <img src="/images/25.png">
    </p></td>
    <td><p align="center">
      <i>Diffuse substrate and glossy coat with an Ashikhmin-Shirley BRDF.</i>
      <img src="/images/24.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>10mm fish-eye lens.</i>
      <img src="/images/21.png">
    </p></td>
    <td><p align="center">
      <i>Subsurface scattering, skin material.</i>
      <img src="/images/19.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>Glass material, lit with an infinite area light backed by a latitude-longitude radiance map.</i>
      <img src="/images/18.png">
    </p></td>
    <td><p align="center">
      <i>Glass material, volumetric path tracing; a bug prevents the environment map from being sampled.</i>
      <img src="/images/17.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>Metal material, lit with an infinite area light backed by a latitude-longitude radiance map.</i>
      <img src="/images/16.png">
    </p></td>
    <td><p align="center">
      <i>Gold material, Torrance-Sparrow microfacets with Trowbridge-Reitz distribution.</i>
      <img src="/images/1.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>Ganesha mesh by Wenzel Jakob; notice the absence of shadows (wrong).</i>
      <img src="/images/2.png">
    </p></td>
    <td><p align="center">
      <i>Participating media. Not working yet.</i>
      <img src="/images/3.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>Cornell box, direct lighting integrator, stratified sampling, 30 samples per pixel.</i>
      <img src="/images/15.png">
    </p></td>
    <td><p align="center">
      <i>Cornell box, path integrator, stratified sampling, 10 samples per pixel.</i>
      <img src="/images/4.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>Plastic material with Torrance-Sparrow microfacet BRDF to simulate glossy hard plastic.</i>
      <img src="/images/5.jpg">
    </p></td>
    <td><p align="center">
      <i>Spherical diffuse area light, matte material with Lambertian reflection.</i>
      <img src="/images/6.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>Spherical diffuse area light, plastic material w/o specular reflection.</i>
      <img src="/images/7.png">
    </p></td>
    <td><p align="center">
      <i>Matte material with Oren-Nayar diffuse reflection.</i>
      <img src="/images/8.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>Matte material with perfect Lambertian diffuse reflection.</i>
      <img src="/images/9.png">
    </p></td>
    <td><p align="center">
      <i>Cover of the PBRT book by Yining Karl Li, perspective camera, 2 point lights, mirror material, specular BRDF.</i>
      <img src="/images/10.jpg">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>Dragon triangle mesh by Christian Schüller, perspective camera, 2 point lights, matte material.</i>
      <img src="/images/11.png">
    </p></td>
    <td><p align="center">
      <i>Cover of the PBRT book by Yining Karl Li, perspective camera, 2 point lights, matte material.</i>
      <img src="/images/12.png">
    </p></td>
  </tr>

  <tr>
    <td><p align="center">
      <i>One triangle of a triangle mesh, matte material.</i>
      <img src="/images/13.png">
    </p></td>
    <td><p align="center">
      <i>First render, implicitly defined sphere.</i>
      <img src="/images/14.png">
    </p></td>
  </tr>
</table>

The sections I've covered are marked next. Each covered section is accompanied by my personal **[notes](notes/)**.

1 Introduction ([notes/1 Introduction](notes/1_Introduction.pdf))
- [x] 1.1 Literate Programming
- [x] 1.2 Photorealistic Rendering and the Ray-Tracing Algorithm
- [x] 1.3 pbrt: System Overview
- [x] 1.4 Parallelization of pbrt
- [x] 1.5 How to Proceed through This Book
- [x] 1.6 Using and Understanding the Code
- [x] 1.7 A Brief History of Physically Based Rendering

2 Geometry and Transformations ([notes/2 Geometry and Transformations](notes/2_Geometry_and_Transformations.pdf))
- [x] 2.1 Coordinate Systems
- [x] 2.2 Vectors
- [x] 2.3 Points
- [x] 2.4 Normals
- [x] 2.5 Rays
- [x] 2.6 Bounding Boxes
- [x] 2.7 Transformations
- [x] 2.8 Applying Transformations
- [ ] 2.9 Animating Transformations
- [x] 2.10 Interactions

3 Shapes ([notes/3 Shapes](notes/3_Shapes.pdf))
- [x] 3.1 Basic Shape Interface
- [x] 3.2 Spheres
- [ ] 3.3 Cylinders
- [ ] 3.4 Disks
- [ ] 3.5 Other Quadrics
- [x] 3.6 Triangle Meshes
- [ ] 3.7 Curves
- [ ] 3.8 Subdivision Surfaces
- [x] 3.9 Managing Rounding Error

4 Primitives and Intersection Acceleration ([notes/4 Primitives and Intersection Acceleration](notes/4_Primitives_and_Intersection_Acceleration.pdf))
- [x] 4.1 Primitive Interface and Geometric Primitives
- [x] 4.2 Aggregates
- [x] 4.3 Bounding Volume Hierarchies
- [ ] 4.4 Kd-Tree Accelerator

5 Color and Radiometry ([notes/5 Color and Radiometry](notes/5_Color_and_Radiometry.pdf))
- [x] 5.1 Spectral Representation
- [x] 5.2 The SampledSpectrum Class
- [x] 5.3 RGBSpectrum Implementation
- [x] 5.4 Radiometry
- [x] 5.5 Working with Radiometric Integrals
- [x] 5.6 Surface Reflection

6 Camera Models ([notes/6 Camera Models](notes/6_Camera_Models.pdf))
- [x] 6.1 Camera Model
- [x] 6.2 Projective Camera Models
- [x] 6.3 Environment Camera
- [x] 6.4 Realistic Cameras

7 Sampling and Reconstruction ([notes/7 Sampling and Reconstruction](notes/7_Sampling_and_Reconstruction.pdf))
- [x] 7.1 Sampling Theory
- [x] 7.2 Sampling Interface
- [x] 7.3 Stratified Sampling
- [ ] 7.4 The Halton Sampler
- [ ] 7.5 (0, 2)-Sequence Sampler
- [ ] 7.6 Maximized Minimal Distance Sampler
- [ ] 7.7 Sobol' Sampler
- [x] 7.8 Image Reconstruction
- [x] 7.9 Film and the Imaging Pipeline

8 Reflection Models ([notes/8 Reflection Models](notes/8_Reflection_Models.pdf))
- [x] 8.1 Basic Interface
- [x] 8.2 Specular Reflection and Transmission
- [x] 8.3 Lambertian Reflection
- [x] 8.4 Microfacet Models
- [ ] 8.5 Fresnel Incidence Effects
- [ ] 8.6 Fourier Basis BSDFs

9 Materials ([notes/9 Materials](notes/9_Materials.pdf))
- [x] 9.1 BSDFs
- [x] 9.2 Material Interface and Implementations
- [ ] 9.3 Bump Mapping

10 Texture ([notes/10 Texture](notes/10_Texture.pdf))
- [x] 10.1 Sampling and Antialiasing
- [x] 10.2 Texture Coordinate Generation
- [x] 10.3 Texture Interface and Basic Textures
- [x] 10.4 Image Texture
- [x] 10.5 Solid and Procedural Texturing
- [ ] 10.6 Noise

11 Volume Scattering
- [x] 11.1 Volume Scattering Processes
- [x] 11.2 Phase Functions
- [x] 11.3 Media
- [x] 11.4 The BSSRDF

12 Light Sources ([notes/12 Light Sources](notes/12_Light_Sources.pdf))
- [x] 12.1 Light Emission
- [x] 12.2 Light Interface
- [x] 12.3 Point Lights
- [ ] 12.4 Distant Lights
- [x] 12.5 Area Lights
- [x] 12.6 Infinite Area Lights

13 Monte Carlo Integration ([notes/13 Monte Carlo Integration](notes/13_Monte_Carlo_Integration.pdf))
- [x] 13.1 Background and Probability Review
- [x] 13.2 The Monte Carlo Estimator
- [x] 13.3 Sampling Random Variables
- [ ] 13.4 Metropolis Sampling
- [x] 13.5 Transforming between Distributions
- [x] 13.6 2D Sampling with Multidimensional Transformations
- [x] 13.7 Russian Roulette and Splitting
- [ ] 13.8 Careful Sample Placement
- [ ] 13.9 Bias
- [x] 13.10 Importance Sampling

14 Light Transport I: Surface Reflection ([notes/14 Light Transport I Surface Reflection](notes/14_Light_Transport_I_Surface_Reflection.pdf))
- [x] 14.1 Sampling Reflection Functions
- [x] 14.2 Sampling Light Sources
- [x] 14.3 Direct Lighting
- [x] 14.4 The Light Transport Equation
- [x] 14.5 Path Tracing

15 Light Transport II: Volume Rendering
- [x] 15.1 The Equation of Transfer
- [x] 15.2 Sampling Volume Scattering
- [x] 15.3 Volumetric Light Transport
- [x] 15.4 Sampling Subsurface Reflection Functions
- [ ] 15.5 Subsurface Scattering Using the Diffusion Equation

16 Light Transport III: Bidirectional Methods
- [ ] 16.1 The Path-Space Measurement Equation
- [ ] 16.2 Stochastic Progressive Photon Mapping
- [ ] 16.3 Bidirectional Path Tracing
- [ ] 16.4 Metropolis Light Transport

17 Retrospective and The Future
- [ ] 17.1 Design Retrospective
- [ ] 17.2 Alternative Hardware Architectures

A Utilities ([notes/Appendix A](notes/Appendix_A.pdf))
- [x] A.1 Main Include File
- [x] A.2 Image File Input and Output
- [x] A.3 Communicating with the User
- [x] A.4 Memory Management
- [x] A.5 Mathematical Routines
- [x] A.6 Parallelism
- [ ] A.7 Statistics

B Scene Description Interface ([notes/Appendix B](notes/Appendix_B.pdf))
- [x] B.1 Parameter Sets
- [x] B.2 Initialization and Rendering Options
- [x] B.3 Scene Definition
- [ ] B.4 Adding New Object Implementations