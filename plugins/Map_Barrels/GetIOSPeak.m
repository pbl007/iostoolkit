% load IOS map
function [mMaskMapT, mMaskMap, mMap] = GetIOSPeak(sFile, nThresh)
%sFile = 'C3_140315.mat';
%nThresh = .00005;

load(sFile)

% Check that ISIdata structure was loaded. If not, then exist
if ~exist('ISIdata', 'var')
    warning(sprintf('GetIOSPeak: %s is not an IOS file.', sFile))
    mMaskMapT = [];
    mMaskMap = [];
    mMap = [];
    return
end

vCLim = [-.0005 .0005];

% invert map (signals are positive)
mMap = double(ISIdata.signalFrame .* -1);

% smooth map
nSigma = 3;
mWin = fspecial('gaussian', nSigma*3, nSigma);
mMap = single(filter2(mWin, mMap, 'same'));

% show map
figure(89)
colormap(hot)
imagesc(mMap)
colorbar
set(gca, 'clim', vCLim)
hTit = title(sFile);
set(hTit, 'interpreter', 'none');
axis equal

% get ROI, which should be drawn conservatively around the region of
% activation
mBW = roipoly;
mMaskMap = mMap .* mBW;

% get values above threshold
mMaskMapT = mMaskMap;
mMaskMapT(mMaskMap <= nThresh) = NaN;

% show masked and thresholded map
surf(mMaskMapT)
shading interp
colorbar
hTit = title(sFile);
set(hTit, 'interpreter', 'none');
axis equal

% save results
%sOutFile = sprintf('%s_MaskedMap.mat', sFile(1:end-4));
%save(sOutFile, 'mMaskMap');
%disp(sprintf('Masked map saved to %s', sOutFile));

return

