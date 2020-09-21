function target = gettarget2DpngNAU(filename,target,im,fig)
%% Inputs
%   filename
%   target
%       xStep_m
%       yStep_m
%       xOffset_m
%       yOffset_m
%       zOffset_m
%       ampAdjust
%       downSample

%% Load the iamge in
tMatrix = imread(filename);
tMatrix = tMatrix(:,:,1);
tMatrix(tMatrix>0) = 1;
tMatrix = ~tMatrix;
tMatrix = fliplr(tMatrix);

%% Crete the image domain
[target.sizeY,target.sizeX] = size(tMatrix);
xAxisT = target.xStep_m * (-(target.sizeX-1)/2 : (target.sizeX-1)/2);
yAxisT = target.yStep_m * (-(target.sizeY-1)/2 : (target.sizeY-1)/2);

xAxisT = xAxisT + target.xOffset_m;
yAxisT = yAxisT + target.yOffset_m;
zAxisT = target.zOffset_m;

[zT,xT,yT] = meshgrid(zAxisT,xAxisT,yAxisT);
target.xyz_m = [xT,yT,zT]; % xPoint x 3 (x-y-z) x yPoint;
target.xyz_m = reshape(permute(target.xyz_m,[1 3 2]),[],3);

indT = rot90(tMatrix,-1)==true;
target.xyz_m = single(target.xyz_m(indT,:));

target.xyz_m = downsample(target.xyz_m,target.downSample);

target.numTarget = size(target.xyz_m,1);
target.amp = ones(target.numTarget,1);

%% Create ideal reflectivity function for png
target.ideal2D = single(zeros(im.numX,im.numY));
for indTarget = 1:target.numTarget
    temp = single(exp(-(target.o_x)^(-2)*(im.x_m-target.xyz_m(indTarget,1)).^2-(target.o_y)^(-2)*(im.y_m-target.xyz_m(indTarget,2)).^2));
    temp = temp/max(temp(:));
    target.ideal2D = target.ideal2D + temp;
end
target.maxpng = max(target.ideal2D(:));
target.ideal2D = target.ideal2D/target.maxpng;

% Adjust the amplitudes of the png to make sure the energy of the png does
% not dwarf the point targets' energy
target.amp = target.amp/target.maxpng*target.ampAdjust;

%% Show the reflectivity function
h = fig.Target2D.h;
mesh(h,im.x_m,im.y_m,target.ideal2D','FaceColor','interp')
view(h,2)
xlabel(h,"x (m)")
ylabel(h,"y (m)")
xlim(h,[im.x_m(1),im.x_m(end)])
ylim(h,[im.y_m(1),im.y_m(end)])
title(h,"Original Reflectivity Function");