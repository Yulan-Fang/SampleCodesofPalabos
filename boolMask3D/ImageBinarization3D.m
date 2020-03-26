clear
clc
%% Clone and padding layers of the image in X axis direction for Palabos 
% 这里的1.png为90×30的图片
% Here the size of 3D.png is 90×30.
Zlength= 30;
%Your intended length of Z axis in Palabos
fprintf('Your intended Z axis length is:%6.2f\n',Zlength)
%% Read image, Binarization, and convert 0 to 1 in the matrix because in Palabos flag 1 for Bouceback nodes.
Read_InputImage=imread('3D.png');  
Read_InputImage=flipud(Read_InputImage);
thresh=graythresh(Read_InputImage);
Binarized_Image=im2bw(Read_InputImage,thresh);
%Binarized_Image=flipud(Binarized_Image);
%Binarized_Image=rot90(B,3);
ConvertValueFlag=~Binarized_Image;
OriginalY=size(ConvertValueFlag,1);
OriginalX=size(ConvertValueFlag,2);
%imshow(ConvertValueFlag);%A preview of 2D domain
fprintf('The thresh of image binarization is:%6.2f\n',thresh)
fprintf('The image size is X=%6.2f, Y=%6.2f\n',OriginalX,OriginalY)
%% Use cat functional in Matlab to create 3D matrix.
MultiCloneofConvertValueFlag=ConvertValueFlag;
for i=1:(Zlength-1) %Padding layers to Multiclone
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
% So we need to make Matlab output with the Z values first.
% The aim of CoordinatorMatrix and Final3Ddata is (Yline,Zcolumn,Xlayers)
% CoordinatorMatrix(Image width,Zlength,Image length)
CoordinatorMatrix=zeros([size(MultiCloneofConvertValueFlag,2),...%Image length
                        size(MultiCloneofConvertValueFlag,3),...%Zlength
                        size(MultiCloneofConvertValueFlag,1)]);%Image width
Final3Ddata=CoordinatorMatrix;%Get the dimension of matrix MulticloneofConvertValueFlag
for YLine=1:size(CoordinatorMatrix,1) %MatrixColumn should be image length, here for example is 270.
    for ZColumn=1:size(CoordinatorMatrix,2)%MatrixLine should be Zlength.
        for ImageLayer=1:size(CoordinatorMatrix,3)%MatrixLayer should be image width.
                  Final3Ddata( YLine,ZColumn,ImageLayer)=...
 MultiCloneofConvertValueFlag(ImageLayer,YLine,ZColumn);
        end
    end
end
FinalY = size(Final3Ddata,1); %Corresponding to X axis in Matlab
FinalZ = size(Final3Ddata,2); %Corresponding to Y axis in Matlab
FinalX = size(Final3Ddata,3);
fprintf('The output 3D matrix size is:\n X=%6.2f, Y=%6.2f(For Palabos Z), Z=%6.2f\n Write the output file in this folder.\n',FinalY,FinalZ,FinalX)
%Write output to 3D.dat with space delimiter.
dlmwrite('3D.dat',Final3Ddata,'delimiter',' ');
isosurface(Final3Ddata)
