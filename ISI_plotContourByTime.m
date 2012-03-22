function h = ISI_plotContourByTime(ISIdata, prmts, cmap)
% h = ISI_plotContourByTime(ISIdata, clim, colormap)
% Plot contour at selected value across all time-stamps in same axes.
%
% Generate temporal map of peak negativity at each pixel (time coded in color)

% Modified by Per M Knutsen , Feb 02 2012
%

if ~isnan(prmts.smoothSigma)
    mWin = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
else
    mWin = NaN;
end

deltaSignal = ISIdata.deltaSignal;

% Smooth all frames
for f = 1:size(ISIdata.deltaSignal, 3)
    if ~isnan(mWin)
        deltaSignal(:,:,f) = single(filter2(mWin, deltaSignal(:,:,f), 'same'));
    end    
end

% Find the minimum value of every pixel and its temporal index
[mMinMap, mMinMapIndx] = min(deltaSignal, [], 3);
[mMaxMap, mMaxMapIndx] = max(deltaSignal, [], 3);

% Remove outliers in mMinMap and mMinMap
nLvl = 6;
mMinMap(mMinMap(:) > (mean(mMinMap(:)) + nLvl*std(mMinMap(:)))) = mean(mMinMap(:));
mMinMap(mMinMap(:) < (mean(mMinMap(:)) - nLvl*std(mMinMap(:)))) = mean(mMinMap(:));

mMaxMap(mMaxMap(:) > (mean(mMaxMap(:)) + nLvl*std(mMaxMap(:)))) = mean(mMaxMap(:));
mMaxMap(mMaxMap(:) < (mean(mMaxMap(:)) - nLvl*std(mMaxMap(:)))) = mean(mMaxMap(:));

% Convert mMinMapIndx and mMaxMapIndx values to seconds
vTime = (((ISIdata.frame_rate/ISIdata.bin_duration):f+((ISIdata.frame_rate/ISIdata.bin_duration)-1)) ...
    / (ISIdata.frame_rate/ISIdata.bin_duration)) ...
    - (ISIdata.nPreStimFrames / (ISIdata.frame_rate/ISIdata.bin_duration)) ...
    - 1;
mMinMapIndx(:) = interp1(1:size(deltaSignal, 3), vTime, mMinMapIndx(:));
mMaxMapIndx(:) = interp1(1:size(deltaSignal, 3), vTime, mMaxMapIndx(:));

% Smooth maps
if ~isnan(mWin)
    mMinMap = single(filter2(mWin, mMinMap, 'same'));
    mMinMapIndx = single(filter2(mWin, mMinMapIndx, 'same'));
    mMaxMap = single(filter2(mWin, mMaxMap, 'same'));
    mMaxMapIndx = single(filter2(mWin, mMaxMapIndx, 'same'));
end

% Plot
hFig = figure;

hAx = subplot(2,2,1);
imagesc(mMinMap)
colorbar
axis image
colormap(cmap)
set(hAx, 'xtick', [], 'ytick', []);%, 'clim', [ISIdata.climAll(1) 0])
title('Minimum intensity projection')

hAx = subplot(2,2,2);
imagesc(mMinMapIndx)
colorbar
axis image
set(hAx, 'xtick', [], 'ytick', [])
colormap(cmap)
title('Minimum time map')

hAx = subplot(2,2,3);
imagesc(mMaxMap)
colorbar
axis image
set(hAx, 'xtick', [], 'ytick', []);%, 'clim', [0 ISIdata.climAll(2)])
colormap(cmap)
title('Maximum intensity projection')

hAx = subplot(2,2,4);
imagesc(mMaxMapIndx)
colorbar
axis image
set(hAx, 'xtick', [], 'ytick', [])
colormap(cmap)
title('Maximum time map')

colormap jet
return
