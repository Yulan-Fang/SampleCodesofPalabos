clear
clc
%% Clone and padding layers of the image
Height= 50;
%Your intended length of height axis in Palabos
fprintf('Your intended height is:%6.2f\n',Height)
%% Read image, Binarization, and convert 0 to 1 in the matrix because in Palabos flag 1 for Bouceback nodes.
Read_InputImage=imread('3D.png');  
Read_InputImage=flipud(Read_InputImage);
thresh=graythresh(Read_InputImage);
Binarized_Image=im2bw(Read_InputImage,thresh);
%Binarized_Image=flipud(Binarized_Image);
%Binarized_Image=rot90(B,3);
ConvertValueFlag=~Binarized_Image;
Image_Width=size(ConvertValueFlag,1);
Image_Length=size(ConvertValueFlag,2);
%imshow(ConvertValueFlag);%A preview of 2D domain
fprintf('The thresh of image binarization is:%6.2f\n',thresh)
fprintf('The image size is length=%6.2f, width=%6.2f\n',Image_Length,Image_Width)
%% Use cat functional in Matlab to create 3D matrix.
MultiCloneofConvertValueFlag=ConvertValueFlag;
for i=1:(Height-1) %Padding layers to Multiclone
    MultiCloneofConvertValueFlag= cat(3,MultiCloneofConvertValueFlag,ConvertValueFlag);
    %Here we get a (Width of picture,Length of picture,Zlength) matrix for Palabos.
end
%% Changing (y,x,Zlength) matrix to (Zlength,y,x) matrix.
% Because in Palabos the data was read and write to domain in a pattern like below:
% for X→Xlength{
%     for Y→Ylength{
%        for Z→Zlength{
%           write to(nx,ny,nz);
% }}}
% So we need to make Matlab output corresponding .
CoordinatorMatrix=zeros([size(MultiCloneofConvertValueFlag,1),...%Image width
                        size(MultiCloneofConvertValueFlag,2),...%Length
                        size(MultiCloneofConvertValueFlag,3)]);%Matlab Matrix Height, which is along Y axis in Palabos
Final3Ddata=CoordinatorMatrix;%Get the dimension of matrix MulticloneofConvertValueFlag
for PLB_Y=1:size(CoordinatorMatrix,3)
    for PLB_Z=1:size(CoordinatorMatrix,2)
        for PLB_X=1:size(CoordinatorMatrix,1)
                  Final3Ddata(PLB_X,PLB_Z,PLB_Y)=...
 MultiCloneofConvertValueFlag(PLB_X,PLB_Z,PLB_Y);
        end
    end
end
Final_row = size(Final3Ddata,1); 
Final_column = size(Final3Ddata,2);
Final_height = size(Final3Ddata,3);
fprintf('The output 3D matrix size is:\n X=%6.2f(For Palabos X), Y=%6.2f(For Palabos Z), Z=%6.2f(For Palabos Y)\nWriting the output file in this folder...\n',Final_row,Final_column,Final_height)
%Write output to 3D.dat with space delimiter.
dlmwrite('3D.dat',Final3Ddata,'delimiter',' ');
isosurface(Final3Ddata)
fprintf('finished\n')
%after showing the figure, we can add descriptions on 3 axis to see, the y
%axis cooresponds to palabos_x, the x axis corresponds to palabos_z, and
%the z axis cooresponds to palabos_y