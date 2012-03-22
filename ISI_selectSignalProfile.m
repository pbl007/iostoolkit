function ISIdata = ISI_selectSignalROI(ISIdata, prmts, cmap)
% ISI_selectSignalProfile   Select a signal profile
%


% Display reference image and let user select ROI for analysis.
mFrame = ISIdata.signalFrame;
mWin = ones(prmts.imgBin);

% Smooth image
if ~isnan(prmts.smoothSigma)
    mWin = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
    mFrame = single(filter2(mWin, mFrame, 'same'));
end

% Read vessel image
sFile = fullfile(prmts.path2dir, prmts.refImage);
mVesselImg = imread(sFile);
mVesselImg = mat2gray(mVesselImg);

% Display vessel image
hFig = figure('Name','Select Signal ROI');

mCols = get(hFig, 'DefaultAxesColorOrder');
mCols = [mCols; mCols ./ 2];

hAx(1) = subplot(1,2,1);
imshow(mVesselImg)
axis image off; hold on
if isfield(ISIdata, 'analysisSignalProfile')
    for i = 1:length(ISIdata.analysisSignalProfile)
        plot(ISIdata.analysisSignalProfile(i).vXi, ISIdata.analysisSignalProfile(i).vYi, 'o-', ...
            'color', mCols(i,:), 'Tag', 'REMOVE')
    end
end

% Display signal image
hAx(2) = subplot(1,2,2);
imagesc(mFrame)
axis image off; hold on
if isfield(ISIdata, 'climAll')
    set(gca, 'Clim', ISIdata.climAll)
end
if isfield(ISIdata, 'analysisSignalProfile')
    for i = 1:length(ISIdata.analysisSignalProfile)
        plot(ISIdata.analysisSignalProfile(i).vXi, ISIdata.analysisSignalProfile(i).vYi, 'o-', ...
            'color', mCols(i,:), 'Tag', 'REMOVE')
    end
end

colormap(cmap)

sAns1 = questdlg('Choose image to draw profile on.', 'Profile Image', ...
    'Vessel Image', 'Average Signal Image', 'Average Signal Image');
if isempty(sAns1), return; end

title('Press ESC to re-use last set profile (in white)')

delete(findobj('Tag', 'REMOVE'))

cXi = {};
cYi = {};
while 1
    if strcmp(sAns1, 'Vessel Image')
        axes(hAx(1))
    else
        axes(hAx(2))
    end
    [cXi{end+1}, cYi{end+1}] = getline;
    axes(hAx(1))
    plot(cXi{end}, cYi{end}, 'o-', 'color', mCols(length(cXi), :))
    axes(hAx(2))
    plot(cXi{end}, cYi{end}, 'o-', 'color', mCols(length(cXi), :))
    sAns = questdlg('Do want to draw one more profile?', 'Profile', ...
        'Yes', 'No', 'No');
    if isempty(sAns) || strcmp(sAns, 'No')
        break
    end
end

% Store results
if ~isempty(cXi)
    ISIdata.analysisSignalProfile = struct('vXi', cXi, 'vYi', cYi);
end

return