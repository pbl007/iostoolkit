% load IOS map
function [mMaskMap, mMap] = GetIOSPeak(sFile, nThresh)

load(sFile) % load .dat file

% Check that ISIdata structure was loaded. If not, then exist
if ~exist('ISIdata', 'var')
    warning(sprintf('GetIOSPeak: %s is not an IOS file.', sFile))
    %mMaskMapT = [];
    mMaskMap = [];
    mMap = [];
    return
end

% invert map (signals are positive)
mMap = double(ISIdata.signalFrame .* -1);

% Get sigma factor for smoothing of maps from GUI
nSigma = str2double(get(findobj(findobj('tag', 'ISIanalysisGUI_fig'), 'tag', 'smooth_sigma'), 'string'));

% Smooth map
if isnumeric(nSigma) && ~isempty(nSigma) && ~isnan(nSigma)
    mWin = fspecial('gaussian', nSigma*3, nSigma);
    mMap = single(filter2(mWin, mMap, 'same'));
end

% Apply manual mask (if set in GUI)
if get(findobj('tag', 'chk_manualmask'), 'value') == 1
    mManualMask = get(findobj('tag', 'btn_setmanualmask'), 'userdata');
    if ~isempty(mManualMask) && all(size(mManualMask.mROI) == size(mMap))
        mMap(~mManualMask.mROI) = NaN;
    end
end

% Determine clim from 1st and 99th percentile of intensity values
%vCLim = [-.0005 .0005];
vCLim = prctile(mMap(:), [.5 99.5]);

% show map
figure(89); clf
colormap(pink)
imagesc(mMap)
colorbar
set(gca, 'clim', vCLim)
hTit = title(sFile);
set(hTit, 'interpreter', 'none');
axis equal off

% get ROI, which should be drawn conservatively around the region of
% activation
mBW = roipoly;
mMaskMap = mMap .* mBW;
mMaskMap(mMaskMap == 0) = NaN;

% get values above threshold
%mMaskMapT = mMaskMap;
%mMaskMapT(mMaskMap <= -nThresh) = NaN;

% show masked and thresholded map
%surf(mMaskMapT)
%shading interp
%colorbar
%hTit = title(sFile);
%set(hTit, 'interpreter', 'none');
%axis equal

return

