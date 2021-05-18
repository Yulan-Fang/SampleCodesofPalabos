% This code is for the slice of the fluid domain
clc;clear;close all
    All=csvread('filename.csv',1);
    vel_mag=All(:,1);
    x=All(:,2);
    y=All(:,3);
    z=All(:,4);
%A velocity Magnitude; B x; C y; D z
%% CSV inherent data
xr = sort(unique(x));
yr = sort(unique(y));
%%
Vel_magnitude=reshape(vel_mag,size(xr,1),size(yr,1));
%%
contourf(rot90(Vel_magnitude(:,:)))
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w')
% filename=['DMD_vel_',num2str(i,'%02.f'),'.png'];
% print(gcf,'-dpng',filename) 
