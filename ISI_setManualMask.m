function tMask = ISI_setManualMask(sFile, hObject)
% M = ISI_SETMANUALMASK(F, H)
% Manually select a mask region on top of vessel image with filepath F. The
% return structure M contains the mask and its vertices. Additionally, mask
% and vertices are stored in the UserData field of object handle H. H can
% optionally be empty.

% Read and display vessel image
mVesselImg = imread(sFile);
mVesselImg = mat2gray(mVesselImg);

% Display vessel image
hFig = figure('Name','Select mask');
hAx = imshow(mVesselImg);
title('Select mask')

[mROImask, vXi, vYi] = roipoly;
imshow(mROImask.*mVesselImg)

% Store results
tMask = struct('mROI', mROImask, 'vXi', vXi, 'vYi', vYi);

if ishandle(hObject)
    set(hObject, 'UserData', tMask)
end

return