## Inspect the stl file
You will need to provide a Makefile to complile and create a tmp folder here.

Things you need to consider:
1) At line 146 you indicate which stl file you want to read.
2) Line 153 and 154 are used to scale the mesh, to convert to LB unit the physical length is scaled by 1/dx
3) You can use paraView to load the vtk file in the tmp folder and the surface_ file just inside the same folder with the cpp file to check.
