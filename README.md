# Description

The repository contains the work testing different schemes in space and in time with the wave equation.
For space discretizations, Finite Difference Method "FDM" (the standard 2nd order and higher space orders) was tested firstly, then Spectral Elements Method "SEM" (based on Finite Element Method) was tested (high polynomial orders in space).
For time discretizations, the LeapFrog scheme and the Modified Equation schemes were tested.

## Organization

The project is organized as follows
* **1D_Tests**: contains the tests of the schemes in time and in space. A first subfolder **Global_time_schemes** contains the files applying the same scheme over the whole mesh. A second subfolder **Mixed_time_schemes** contains the tests with multiple time schemes for a mesh containing one tiny element in comparison to the rest.
  
* **2D_Tests**: contains the tests in 2D with SEM in space and modified equation high orders in time. It also contains the codes to transform the mesh from gmsh ".msh" files to readable data in matlab ad creates the corresponding mass and stiffness matrices M and K.
  
* **Archive**: contains the different tests made during work, with various intermediate results with the "unpolished" versions of the codes in **1D_Tests** and **2D_Tests** 
