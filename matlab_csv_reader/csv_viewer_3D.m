clc;clear;close all
file_read=dir('*.csv');
filename={file_read.name};
file_length=length(file_read);

i=91 % select the ith file to read
    filename1=file_read(i).name;
    All=csvread(filename1,1);
    velocitynorm=All(:,1);
    u=All(:,2);
    v=All(:,3);
    w=All(:,4);
    vor_u=All(:,5);
    vor_v=All(:,6);
    vor_w=All(:,7);
    vor_magnitude=All(:,8);

%A velocitynorm; B u; C v; D w; E vor_u; F vor_v; G vor_w; H vor_magnitude
%% CSV inherent data
load xdata.mat;% I already save the x,y columns as mat file in the same directory
load ydata.mat;
xr = sort(unique(x));
yr = sort(unique(y));

Velocity_x=reshape(u,size(xr,1),size(yr,1),size(yr,1));

% I can reshape the column like this because the y and z length of my simulation domain are the same.

Velocity_y=reshape(v,size(xr,1),size(yr,1),size(yr,1));
Velocity_z=reshape(w,size(xr,1),size(yr,1),size(yr,1));
Velocitynorm=reshape(w,size(xr,1),size(yr,1),size(yr,1));
Vor_u=reshape(vor_u,size(xr,1),size(yr,1),size(yr,1));
Vor_v=reshape(vor_v,size(xr,1),size(yr,1),size(yr,1));
Vor_w=reshape(vor_w,size(xr,1),size(yr,1),size(yr,1));
Vor_magnitude=reshape(vor_magnitude,size(xr,1),size(yr,1),size(yr,1));

%%
vel_magnitude=sqrt(u.^2+v.^2+w.^2);
Vel_magnitude=reshape(vel_magnitude,size(xr,1),size(yr,1),size(yr,1));
%%
contourf(flipud(rot90(Vel_magnitude(:,:,66))))
%imagesc(flipud(rot90(Vel_magnitude(:,:,66))))
colormapeditor
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
% filename=['FromCsv_vel_',num2str(i,'%02.f'),'.png'];
% print(gcf,'-dpng',filename) 
