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

nScale = 1;

mVesselImg = imread(sVesselFile);
mVesselImg = imresize(mVesselImg, nScale, 'bicubic');
mVesselImgOrig = mVesselImg;
mVesselImg = mat2gray(mVesselImg); %don't really need to normalize, but do it for clarity

% smooth
nSigma = str2double(get(findobj(findobj('tag', 'ISIanalysisGUI_fig'), 'tag', 'smooth_sigma'), 'string'));
if isnumeric(nSigma) && ~isempty(nSigma) && ~isnan(nSigma)
    mWin = fspecial('gaussian', nSigma*3, nSigma);
    mVesselImg = single(filter2(mWin, mVesselImg, 'same'));
end


% subtract a centered gaussian
mWin = gausswin(size(mVesselImg, 1)) * gausswin(size(mVesselImg,2))';

mVesselImg = mVesselImg - mWin./(4./median(mVesselImg(:)));

% remove zeros and ones
mVesselImg(mVesselImg == 0) = NaN;
mVesselImg(mVesselImg == 1) = NaN;

% drop 5% percentile, both ends
vMinMax = prctile(mVesselImg(:), [1 99]);

% normalize image
%mVesselImg(mVesselImg <= vMinMax(1)) = NaN;
%mVesselImg(mVesselImg >= vMinMax(2)) = NaN;
%mVesselImg = mVesselImg - vMinMax(1);
%mVesselImg = mVesselImg ./ max(mVesselImg(:));


% Apply manual mask if it exists (check GUI)
hBtn = findobj('Tag', 'btn_setmanualmask'); % handle to 'Set' button
if ishandle(hBtn)
    tMask = get(hBtn, 'UserData');
    if ~isempty(tMask)
        tMask.mROI = imresize(tMask.mROI, nScale, 'nearest');
        mVesselImg(~tMask.mROI) = NaN;
    end
end


vmask = im2bw(mVesselImg, vthresh);

vmaskrgb = repmat(vmask, [1 1 3]);
vmaskrgb(:,:,2:3) = 0;

% TODO:
% Single white dots surrounded by only black dots become black
% Single black dots surrounded by only black dots become white


hFig = findobj('tag', 'ISI_maskPreview');
if isempty(hFig)
    hFig = figure;
    set(hFig, 'tag', 'ISI_maskPreview')
else
    figure(hFig)
end

subplot(1,3,1)
imagesc(mVesselImgOrig);
axis image off;
title('Original')

subplot(1,3,2)
imagesc(vmask);
axis image off;
title(['Mask threshold = ' num2str(thresh)]);

subplot(1,3,3)
hv = imshow(mVesselImg);  
axis image; hold on;
hmask = imshow(vmask);
set(hmask,'alphadata', 0.3, 'CData', vmaskrgb);
title('Overlay')
hold off

return