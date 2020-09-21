function target = gettarget2DrandNAU(target,im,fig)
%% Inputs
%   filename
%   target
%       numTargetMax
%       zOffset_m
%       o_x
%       o_y
%       ampMin
%       ampMax

%% Create the target locations and amplitudes
target.numTarget = randi(target.numTargetMax);
target.xyz_m = single([im.x_m(1) + (im.x_m(end)-im.x_m(1))*rand(target.numTarget,1),im.x_m(1) + (im.y_m(end)-im.y_m(1))*rand(target.numTarget,1),target.zOffset_m*ones(target.numTarget,1)]);
target.amp = target.ampMin + (target.ampMax-target.ampMin)*rand(1,target.numTarget);

fail = true;
tic;

while fail
    R_same = pdist2(target.xyz_m,target.xyz_m) + 1e3*eye(target.numTarget);
    indGood = min(R_same,[],1)>(target.o_x*2);
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

%% Show the reflectivity function
h = fig.Target2D.h;
mesh(h,im.x_m,im.y_m,target.ideal2D','FaceColor','interp')
view(h,2)
xlabel(h,"x (m)")
ylabel(h,"y (m)")
xlim(h,[im.x_m(1),im.x_m(end)])
ylim(h,[im.y_m(1),im.y_m(end)])
title(h,"Original Reflectivity Function, " + target.numTarget + " targets");
