function ISI_viewMask(sFilename, thresh)
% View the mask on top of the original vessel file
% Ssed to compare actual coverage, see what areas are being masked out.

if nargin==0
    sVesselFile = fullfile(basepath, subpath,sFilename);
    thresh = 0.45;
else
    sVesselFile = sFilename;
end
 
vthresh = thresh;

mVesselImg = imread(sVesselFile);

% Apply manual mask if it exists (check GUI)
hBtn = findobj('Tag', 'btn_setmanualmask'); % handle to 'Set' button
if ishandle(hBtn)
    tMask = get(hBtn, 'UserData');
    if ~isempty(tMask)
        mVesselImg(~tMask.mROI) = NaN;
    end
end

mVesselImg = mat2gray(mVesselImg); %don't really need to normalize, but do it for clarity
vmask = im2bw(mVesselImg, vthresh);
vmaskrgb = repmat(vmask, [1 1 3]);
vmaskrgb(:,:,2:3) = 0;

figure
hv = imshow(mVesselImg);  
axis image; hold on;
hmask = imshow(vmask);
set(hmask,'alphadata', 0.3, 'CData', vmaskrgb);

title(['mask threshold = ' num2str(thresh)]);

return