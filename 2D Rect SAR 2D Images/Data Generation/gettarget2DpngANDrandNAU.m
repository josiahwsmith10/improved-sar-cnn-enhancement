function target = gettarget2DpngANDrandNAU(filename,target,im,fig)
%% Inputs
%   filename
%   target
%       numTargetMax
%       xStep_m
%       yStep_m
%       xOffset_m
%       yOffset_m
%       zOffset_m
%       o_x
%       o_y
%       ampMin
%       ampMax
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
png.xyz_m = [xT,yT,zT]; % xPoint x 3 (x-y-z) x yPoint;
png.xyz_m = reshape(permute(png.xyz_m,[1 3 2]),[],3);

indT = rot90(tMatrix,-1)==true;
png.xyz_m = single(png.xyz_m(indT,:));

png.xyz_m = downsample(png.xyz_m,target.downSample);

png.numTarget = size(png.xyz_m,1);
png.amp = ones(png.numTarget,1);

%% Create ideal reflectivity function for png
png.ideal2D = single(zeros(im.numX,im.numY));
for indTarget = 1:png.numTarget
    temp = single(exp(-(target.o_x)^(-2)*(im.x_m-png.xyz_m(indTarget,1)).^2-(target.o_y)^(-2)*(im.y_m-png.xyz_m(indTarget,2)).^2));
    temp = temp/max(temp(:));
    png.ideal2D = png.ideal2D + temp;
end
png.maxpng = max(png.ideal2D(:));
png.ideal2D = png.ideal2D/png.maxpng;

% Adjust the amplitudes of the png to make sure the energy of the png does
% not dwarf the point targets' energy
png.amp = png.amp/png.maxpng*target.ampAdjust;

%% Create the target locations and amplitudes
target.numTarget = randi(target.numTargetMax);
fail = true;

tic
target.xyz_m = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget,1),im.x_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget,1),target.zOffset_m*ones(target.numTarget,1)]);
while fail
    target.amp = target.ampMin + (target.ampMax-target.ampMin)*rand(target.numTarget,1);
    R = pdist2(png.xyz_m,target.xyz_m);
    R_same = pdist2(target.xyz_m,target.xyz_m) + 1e3*eye(target.numTarget);
    indGood = min(R,[],1)>(target.o_x*2) & min(R_same,[],1)>(target.o_x*2);
    numGood = sum(indGood);
    if numGood == target.numTarget
        fail = false;
    else
        xyz_m_good = target.xyz_m(indGood,:);
        xyz_m_new = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget-size(xyz_m_good,1),1),im.x_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget-size(xyz_m_good,1),1),target.zOffset_m*ones(target.numTarget-size(xyz_m_good,1),1)]);
        target.xyz_m = cat(1,xyz_m_good,xyz_m_new);
        
        if toc > 10
            indGood = min(R,[],1)>(target.o_x*2);
            target.xyz_m = target.xyz_m(indGood,:);
            target.amp = target.amp(indGood);
            target.numTarget = sum(indGood);
            warning("Could not place targets correctly in 10s, reducing number of targets")
            break;
        end
    end
end

%% Create the ideal reflectivity function
target.ideal2D = single(zeros(im.numX,im.numY));
for indTarget = 1:target.numTarget
    temp = single(exp(-(target.o_x)^(-2)*(im.x_m-target.xyz_m(indTarget,1)).^2-(target.o_y)^(-2)*(im.y_m-target.xyz_m(indTarget,2)).^2));
    temp = temp*target.amp(indTarget)/max(temp(:));
    target.ideal2D = target.ideal2D + temp;
end
if abs(max(target.ideal2D(:))-max(target.amp)) > 3e-2
    warning("Amplitude error! The ideal image is distorted!")
end
target.ideal2D(target.ideal2D>1) = 1;
target.ideal2D(target.ideal2D<0) = 0;

target.ideal2D = target.ideal2D + png.ideal2D;
target.xyz_m = cat(1,target.xyz_m,png.xyz_m);
target.amp = cat(1,target.amp,png.amp)';
target.numTarget = target.numTarget + png.numTarget;

%% Show the reflectivity function
h = fig.Target2D.h;
mesh(h,im.x_m,im.y_m,target.ideal2D','FaceColor','interp')
view(h,2)
xlabel(h,"x (m)")
ylabel(h,"y (m)")
xlim(h,[im.x_m(1),im.x_m(end)])
ylim(h,[im.y_m(1),im.y_m(end)])
title(h,"Original Reflectivity Function, " + target.numTarget + " targets");