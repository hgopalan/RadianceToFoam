This utlity maps the radiance output to OpenFOAM for thermal simulations. The radiance output has to be transformed similar to setting up the boundary condition for timeVaryingMappedFixedValue. Works for both steady and unsteady cases. Developed for OpenFoam 2.4.X . May need modifications for newer versions. Long-wave is not currently included. However, it is straight-forward to modify the code and read long-wave information from object registry. 

**** The utility is experimental and crashes can be expected ***** 

- Venugopalan and Harish 
