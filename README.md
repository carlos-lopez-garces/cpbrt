This repository is the result of studying Physically Based Rendering: From Theory to Implementation, by Pharr, Jakob, and Humphreys ([online edition](https://pbrt.org/)).

The 1st objective is to study the sections that lead to a minimal path tracer and write the corresponding code.

The sections I've covered are marked next. Each covered section is accompanied by my personal **[notes](notes/)**.

1 Introduction ([notes/1 Introduction](notes/1_Introduction.pdf))
- [x] 1.1 Literate Programming
- [x] 1.2 Photorealistic Rendering and the Ray-Tracing Algorithm
- [x] 1.3 pbrt: System Overview
- [ ] 1.4 Parallelization of pbrt
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
- [ ] 2.10 Interactions

3 Shapes ([notes/3 Shapes](notes/3_Shapes.pdf))
- [x] 3.1 Basic Shape Interface
- [ ] 3.2 Spheres
- [ ] 3.3 Cylinders
- [ ] 3.4 Disks
- [ ] 3.5 Other Quadrics
- [ ] 3.6 Triangle Meshes
- [ ] 3.7 Curves
- [ ] 3.8 Subdivision Surfaces
- [ ] 3.9 Managing Rounding Error

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
- [ ] 6.3 Environment Camera
- [ ] 6.4 Realistic Cameras

7 Sampling and Reconstruction ([notes/7 Sampling and Reconstruction](notes/7_Sampling_and_Reconstruction.pdf))
- [x] 7.1 Sampling Theory
- [x] 7.2 Sampling Interface
- [x] 7.3 Stratified Sampling
- [ ] 7.4 The Halton Sampler
- [ ] 7.5 (0, 2)-Sequence Sampler
- [ ] 7.6 Maximized Minimal Distance Sampler
- [ ] 7.7 Sobol’ Sampler
- [x] 7.8 Image Reconstruction
- [x] 7.9 Film and the Imaging Pipeline

8 Reflection Models ([notes/8 Reflection Models](notes/8_Reflection_Models.pdf))
- [x] 8.1 Basic Interface
- [ ] 8.2 Specular Reflection and Transmission
- [x] 8.3 Lambertian Reflection
- [ ] 8.4 Microfacet Models
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
- [ ] 10.4 Image Texture
- [ ] 10.5 Solid and Procedural Texturing
- [ ] 10.6 Noise

11 Volume Scattering
- [ ] 11.1 Volume Scattering Processes
- [ ] 11.2 Phase Functions
- [ ] 11.3 Media
- [ ] 11.4 The BSSRDF

12 Light Sources ([notes/12 Light Sources](notes/12_Light_Sources.pdf))
- [x] 12.1 Light Emission
- [x] 12.2 Light Interface
- [x] 12.3 Point Lights
- [ ] 12.4 Distant Lights
- [ ] 12.5 Area Lights
- [ ] 12.6 Infinite Area Lights

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
- [ ] 15.1 The Equation of Transfer
- [ ] 15.2 Sampling Volume Scattering
- [ ] 15.3 Volumetric Light Transport
- [ ] 15.4 Sampling Subsurface Reflection Functions
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
- [ ] A.2 Image File Input and Output
- [ ] A.3 Communicating with the User
- [x] A.4 Memory Management
- [x] A.5 Mathematical Routines
- [ ] A.6 Parallelism
- [ ] A.7 Statistics

B Scene Description Interface ([notes/Appendix B](notes/Appendix_B.pdf))
- [x] B.1 Parameter Sets
- [ ] B.2 Initialization and Rendering Options
- [x] B.3 Scene Definition
- [ ] B.4 Adding New Object Implementations