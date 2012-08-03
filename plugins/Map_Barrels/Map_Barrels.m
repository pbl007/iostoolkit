function ISI_MapBarrels(varargin)
% Compute combined barrels maps from IOS signal frames
%
% Graphical outputs:
%   - a combined barrels map
%   - uniquely identified barrel regions
%   - barrel map representations overlaid on vessel image
%   - map of normalized distance within barrels from barrel edges
%
% Saved outputs:
%   *Files are saved to the same directory where .mat files reside.*
%
%   tMaps.mat       Contains a structure with fields;
%                       .id
%                       .raw_map
%                       .thresholded_map
%
%   IOSPeaks.mat    Contains a structure with interactively collected data
%
% Pre-requisites:
% As a pre-requisite to running the plugin, you will first need to save
% signal frames to .mat files through the main window. This can be done as
% a batch job on the entire directory.
%
% Before running this function, the IOS_analysisGUI must already be running
% and loaded with relevant parameters. At the very minimum, the Directory
% and Vessel File fields in the GUI must not be empty.
%
% Notes:
% The plugin saves its results once the all-barrel map has been computed.
% If a saved map is found on disk (in the same directory where .mat files
% are located), then the map-generation part will be skipped and only
% map-representations shown.
% 
% By default, all .mat files on disk are used. If some files are to be
% excluded (e.g. if a barrel was imaged several times), an exclusion list
% can be specified in the file ExcludeFiles.m (which needs to reside in the
% same folder as .mat files). See code below for format of this file.
% 
% Created by;
% Per M Knutsen <pmknutsen@gmail.com>, June 2012.
%

handles = varargin{3};

% Fixed variables
nThresh = 0.00005;

nRows = 2;
nCols = 3;

% Get path from GUI
sPath = get(findobj(findobj('tag', 'ISIanalysisGUI_fig'), 'tag', 'path'), 'string');

% Get vessel image path and filename from GUI
sVesselImage = fullfile(sPath, get(findobj(findobj('tag', 'ISIanalysisGUI_fig'), 'tag', 'vessel_filename'), 'string'));

% Get sigma factor for smoothing of maps from GUI
nSigma = str2double(get(findobj(findobj('tag', 'ISIanalysisGUI_fig'), 'tag', 'smooth_sigma'), 'string'));

% Load the file exclusion list
% This file is be default called ExcludeFiles.m and contains a cell
% called cExclude (as Matlab code). The cell contains filenames that will NOT be included in
% the all-barrel map. Example;
%    cExclude = {    'C1_173543', ...
%                    'C1_174151', ...
%                    'C4_175521' };
sExcludeFile = fullfile(sPath, 'ExcludeFiles.m');
if exist(sExcludeFile, 'file')
    sPwd = pwd;
    cd(sPath)
    ExcludeFiles
    cd(sPwd)
    % If cExclude variable was not loaded, create an empty cell to prevent
    % errors further down
    if ~exist('cExclude', 'var')
        cExclude = {};
    end
else
    cExclude = {};
end

% Check if IOSPeaks.mat file exists (i.e. has already been generated). If
% to, skip loading of .mat files and go straight to the graphics generation
% part
if exist(fullfile(sPath, 'IOSPeaks.mat'), 'file')
    disp(sprintf('** GetIOSPeaksWrapper says:\n** IOSPeaks.mat will be loaded from disk.\n** To generate a new IOSPeaks.mat file, delete the file:\n** %s', fullfile(sPath, 'IOSPeaks.mat')));
    load(fullfile(sPath, 'IOSPeaks.mat'))
    sList = tIOSPeaks(1).sList;
    cMaskMapT = tIOSPeaks(1).cMaskMapT;
    cMaskMap = tIOSPeaks(1).cMaskMap;
    cMap = tIOSPeaks(1).cMap;
    
else
    % If IOSPeaks.mat does not exist, we will generate it by iterating
    % across all .mat files and interactively draw a boundary/ROI around
    % each putative barrel manually
    
    % Get list of all .mat files (mapped whiskers)
    sPwd = pwd;
    cd(sPath)
    sList = dir('*.mat');

    % get all masked maps (interactively since we need to draw the ROI for each file)
    cMaskMapT = cell({});
    cMaskMap = cell({});
    cMap = cell({});
    bHaveData = false;
    for i = 1:length(sList)
        [mMaskMapT, mMaskMap, mMap] = GetIOSPeak(sList(i).name, nThresh);
        if isempty(mMaskMapT)
            continue
        end
        cMaskMapT{end+1} = mMaskMapT;
        cMaskMap{end+1} = mMaskMap;
        cMap{end+1} = mMap;
        if ~exist('sListUse')
            sListUse = sList(i);
        else
            sListUse(end+1) = sList(i);
        end
        bHaveData = true;
    end

    if ~bHaveData
        error('No data found in current IOS folder')
    end
    
    % generate structure to hold data
    tIOSPeaks = struct([]);
    tIOSPeaks(1).cMaskMap = cMaskMap;
    tIOSPeaks(1).cMaskMapT = cMaskMapT;
    tIOSPeaks(1).cMap = cMap;
    tIOSPeaks(1).sList = sListUse;

    % Save IOSPeaks.mat file
    save IOSPeaks.mat tIOSPeaks
    
    cd(sPwd)
end

% Remove empty indices in cMaskMapT cMaskMap cMap
vIndxRem = [];
for i = 1:length(cMap)
    if isempty(cMap{i}) && isempty(cMaskMapT{i}) && isempty(cMaskMap{i})
        vIndxRem(end+1) = i;
    end
end
cMap(vIndxRem) = [];
cMaskMapT(vIndxRem) = [];
cMaskMap(vIndxRem) = [];

% Generate a structure with all individual maps
%   .id
%   .raw_map
%   .thresholded_map
tMaps = struct('id', [], 'raw_map', [], 'thresholded_map', [], 'thresh', []);
for i = 1:size(cMaskMapT, 2)
    tMaps(i).id = sList(i).name(1:2);
    tMaps(i).raw_map = cMap{i};
    tMaps(i).thresholded_map = cMaskMapT{i};
    tMaps(i).thresh = nThresh;
end

% Generate the max projection all-barrels mapcMaskMap
mMaskMapAll = zeros(size(cMaskMapT{1}));
mIdentityMap = zeros(size(cMaskMapT{1}));
vY = []; vX = []; cTxt = cell({});
for i = 1:length(cMaskMapT)
    if any(strcmp(sList(i).name(1:end-4), cExclude))
        continue
    end
    mMap = cMaskMapT{i};
    %mMap = cMap{i};
    mMap(isnan(mMap)) = 0;
    
    % smooth map
    if isnumeric(nSigma) && ~isempty(nSigma) && ~isnan(nSigma)
        mWin = fspecial('gaussian', nSigma*3, nSigma);
        mMap = single(filter2(mWin, mMap, 'same'));
    end

    % normalize
    vUniq = sort(unique(mMap(:)));
    if length(vUniq) == 1, continue, end
    mMap = mMap - vUniq(2);
    mMap = mMap ./ max(max(mMap));
    if all(isnan(mMap(:))), continue, end
    
    % Drop values below 10th percentile (to get rid of outliers)
    mMap(mMap(:) < .3) = 0;

    mMaskMapAll = max(mMap, mMaskMapAll);
    [nMax, nMaxInd] = max(mMap(:));

    % Find indices in mMaskMapAll that were used from current map
    % Values that are ambiguous are assigned NaNs
    mIdentityMap(mMaskMapAll == mMap) = i+1;
    
    % info for adding text labels later
    [nX, nY] = ind2sub(size(mMap), nMaxInd);
    vX(end+1) = nX;
    vY(end+1) = nY;
    cTxt{end+1} = sList(i).name(1:2);

end

% Add combined maps to tMaps
tMaps(end+1).id = 'All_CombinedMap';
tMaps(end).raw_map = [];
tMaps(end).thresholded_map = mMaskMapAll;
tMaps(end).thresh = NaN;
tMaps(end+1).id = 'All_IdentityMap';
tMaps(end).raw_map = [];
tMaps(end).thresholded_map = mIdentityMap;
tMaps(end).thresh = NaN;

% Save tMaps and then dump it
sPwd = pwd;
cd(sPath)
save tMaps.mat tMaps
cd(sPwd)
clear tMaps

mMaskMapAll(mMaskMapAll(:) == 0) = -.5;
mMaskMapAll(mMaskMapAll(:) == NaN) = -.5;

% Display all-barrel map (filtered)
hMainFig = figure;
nSubNo = 1;
subplot(nRows,nCols,nSubNo)
imagesc(mMaskMapAll .* -1)
colormap gray
axis image
title('Combined Barrel Map')
vCLim = [-1 -.5];
set(gca,'clim', vCLim)
xlabel('px'); ylabel('px')

% Text labels
for i = 1:length(vX)
    hTxt = text(vY(i), vX(i), cTxt{i});
    set(hTxt, 'color', 'w', 'horizontalAlignment','center', 'verticalAlignment', 'middle');
end

% Generate mask for mIdentityMap
mMask = ones(size(mIdentityMap));
mMask( (mMaskMapAll .* -1) < vCLim(1) | (mMaskMapAll .* -1) > vCLim(2) ) = 0;

% Display identity map
mIdentityMap = mIdentityMap .* mMask;
figure(hMainFig)
nSubNo = nSubNo + 1;
subplot(nRows,nCols,nSubNo)
mMap = jet(max(mIdentityMap(:)));
mMap(1,:) = [1 1 1]; % white background
subimage(mIdentityMap, mMap)
axis image
title('Barrel Regions w/Centroids (peaks)')
xlabel('px'); ylabel('px')

% Compute centroid of each whisker representation and plot on top of
% identify map
hold on
mCentroids = [];
for i = 1:length(cMaskMapT)
    vIndx = find(mIdentityMap(:) == i);
    [vI, vJ] = ind2sub(size(mIdentityMap), vIndx);
    mCentroids(i, :) = [median(vI) nanmedian(vJ)];
end
scatter(mCentroids(:, 2), mCentroids(:, 1), 'k+')

% Plot barrel map on vessel image
figure(hMainFig)
nSubNo = nSubNo + 1;
subplot(nRows,nCols,nSubNo)
mBlue = imread(sVesselImage);
mBlue = double(mBlue(:,:,1));
mBlue = mBlue - min(mBlue(:));

% Mask vessel image
if get(handles.chk_manualmask, 'value')
    tMask = get(handles.btn_setmanualmask, 'userdata');
    if ~isempty(tMask)
        mBlue(~tMask.mROI) = max(mBlue(:));
    end
end

mBlue = uint8 ( (mBlue ./ max(mBlue(:))) .* 255 );
mBlue = repmat(mBlue, [1 1 3]);

h1 = imshow(mBlue); % uint8 [H W 3]    0 - 255
hold on

mNorm = double(mMaskMapAll) .* -1;

mNorm(mNorm(:) < vCLim(1)) = vCLim(1);
mNorm(mNorm(:) > vCLim(2)) = vCLim(2);

mNorm = mNorm - min(mNorm(:));
mNorm = mNorm ./ max(mNorm(:));
mMap = jet(255)
h2 = subimage(mNorm.*255, mMap); % double [H W]                  0.0001 - 1
alphamap = zeros(size(mNorm));
alphamap(mNorm < .8) =.3;%  on dual-screen setups on Linux, Matlab crashed when transparency is used
set(h2, 'alphadata', alphamap, 'alphaDataMapping','none');
title('Barrel Map overlay on vessel image')

% Plot identity map and centroids on vessel images
figure(hMainFig)
nSubNo = nSubNo + 1;
subplot(nRows,nCols,nSubNo)
h1 = imshow(mBlue); % uint8 [H W 3]    0 - 255
hold on

mNorm = mIdentityMap;
mNorm = mNorm - min(mNorm(:));
mNorm = mNorm ./ max(mNorm(:));

mMap = jet(255);
h2 = subimage(mNorm.*255, mMap); % double [H W]                  0.0001 - 1
alphamap = zeros(size(mNorm));
alphamap(mNorm > 0) = .4; % .4
set(h2, 'alphadata', alphamap, 'alphaDataMapping','none');
scatter(mCentroids(:, 2), mCentroids(:, 1), 'k+')
title('Barrel Regions overlay on vessel image')

% Get contours of regions
figure(hMainFig)
nSubNo = nSubNo + 1;
hAx = subplot(nRows,nCols,nSubNo);
hFigTemp = figure;
mColMap = jet(length(cMaskMap));
mColMap(1,:) = [1 1 1];
mDistanceMap = zeros(size(mIdentityMap)); % map encoding distance from edge of barrels
for i = 1:length(cMaskMap)
    mTempMap = zeros(size(mIdentityMap));
    mTempMap(mIdentityMap == i) = 1;
    
    % Get outline of region with contour() function
    figure(hFigTemp)
    [mC, vH] = contour(mTempMap, 1);

    % Fill region with values ranging from 0 to 1, where;
    %   0 is the edge
    %   1 is the center
    %   intermediate values denote distance from nearest edge
    
    % for all internal pixels, compute eucledian distances to all boundary
    % pixels
    [vI, vJ] = find(mTempMap); % all internal (and boundary) pixels
    if isempty(vI), continue, end
    vMinDists = zeros(size(vI));
    for j = 1:length(vI)
        vMinDists(j) = min( sqrt( [mC(1,:) - vJ(j)].^2 + [mC(2,:) - vI(j)].^2 ) );
    end
    % Normalize distance vector (0.01 - 1)
    vMinDists = [vMinDists - min(vMinDists)] ./ max(vMinDists - min(vMinDists));
    % Replace distance values in temporary matrix
    mTempMap(mTempMap == 1) = vMinDists;
    % Merge single distance map with all other barrels
    mDistanceMap = mDistanceMap + mTempMap;
    
    % Plot into contours figure
    axes(hAx)
    hold on
    plot(mC(1,2:end), mC(2,2:end), 'color', mColMap(i, :))
    
end
close(hFigTemp)
axes(hAx)
axis equal; box on
set(gca, 'xlim', [0 size(mTempMap, 1)], 'ylim', [0 size(mTempMap, 2)], 'ydir', 'reverse')
xlabel('px'); ylabel('px')
title('Barrel Region Contours')

% Show distance map
figure(hMainFig);
nSubNo = nSubNo + 1;
subplot(nRows,nCols,nSubNo)
imagesc(mDistanceMap)
set(hMainFig, 'userdata', mDistanceMap)
axis equal; box on
set(gca, 'xlim', [0 size(mTempMap, 1)], 'ylim', [0 size(mTempMap, 2)])
xlabel('px'); ylabel('px')
%colormap pink
hColAx = colorbar('east');
ylabel(hColAx, 'Distance from edge')
title('Distance from barrel edge (normalized)')
hold on
[mC, hCont] = contour(mDistanceMap,.25,'g');
set(hCont,'tag','cnt')
% Create a slider in figure to set contour level
% slider removes current contours and shows new with new values
uicontrol('Style', 'slider', 'units', 'norm', 'position', [.05 0.01 .9, .03], ...
    'min', 0, 'max', 1, ...
    'value', .25, ...
    'callback', 'h=findobj(''tag'',''cnt'');delete(h);m=get(gcf,''userdata'');[c,h]=contour(m,get(gco,''value''),''g'');set(h,''tag'',''cnt'');');
mColMap = pink;
%mColMap(1,:) = [1 1 1];
colormap(mColMap)

header(sprintf('%s', sPath), 14)

return




function y = header(x, fontsize)
%HEADER Puts a text string on top of the current figure
if nargin < 2; fontsize = 20; end;
t = findobj(gcf, 'Tag', 'header');
if length(t) == 0;
	ax = gca;
	axes('Position', [0 0 1 1], 'Visible', 'off');
	t = text(0.5, 0.98, x, ...
				'Units', 'normalized', ...
				'VerticalAlignment', 'top', ...
				'HorizontalAlignment', 'center', ...
				'Tag', 'header', ...
				'FontSize', fontsize);
	axes(ax);
else;
	set(t, 'String', x, 'FontSize', fontsize);
end;
if nargout > 0; y = t; end;
return



