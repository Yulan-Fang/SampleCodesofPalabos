Normally, Paraview can save the simulation data of Palabos as *.csv files, we may save the slice or the entire domain's data. In those csv files, the data is listed as columns and Matlab can easily handle them.
#### csv_viewer_2D
This code reads slice data and reshape the data from x and y dimensions, then use contourf to visualize it.
