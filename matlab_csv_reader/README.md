Normally, Paraview can save the simulation data of Palabos as *.csv files, we may save the slice or the entire domain's data. In those csv files, the data is listed as columns and Matlab can easily handle them.
#### csv_viewer_2D
This code reads slice data and reshape the data from x and y dimensions, then use contourf to visualize it.
#### csv_viewer_3D
Because it reads 3D data, so the reshape function should be including coordinates of z axis, in this example file I can use y direction twice because my simulation domain's y and z length are the same. If you want the code automatically reads all files, you can simply modify the 'i' to be for 'i=1:file_length' and add the 'end' in the end of the code.
