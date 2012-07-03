function ISIdata = ISI_selectSignalROI(ISIdata, prmts, cmap)
% ISI_selectSignalROI
% Manually select regions of interest on top of vessel or signal images.

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
vPos = get(hFig, 'position');
set(hFig, 'position', [vPos(1:2) 700 350]);

mCols = get(hFig, 'DefaultAxesColorOrder');
mCols = [mCols; mCols ./ 2];

hAx(1) = subplot(1,2,1);
imshow(mVesselImg)
axis image off; hold on
if isfield(ISIdata, 'analysisSignalROI')
    for i = 1:length(ISIdata.analysisSignalROI)
        plot(ISIdata.analysisSignalROI(i).vXi, ISIdata.analysisSignalROI(i).vYi, 'o-', ...
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
if isfield(ISIdata, 'analysisSignalROI')
    for i = 1:length(ISIdata.analysisSignalROI)
        plot(ISIdata.analysisSignalROI(i).vXi, ISIdata.analysisSignalROI(i).vYi, 'o-', ...
            'color', mCols(i,:), 'Tag', 'REMOVE')
    end
end

colormap(cmap)

sAns = questdlg('Choose image to draw ROI on', 'ROI Image', ...
    'Vessel Image', 'Average Signal Image', 'Average Signal Image');
if isempty(sAns), return; end

title('Press ESC to re-use last set ROI (in white)')

delete(findobj('Tag', 'REMOVE'))

cROImask = {};
cXi = {};
cYi = {};
while 1
    if strcmp(sAns, 'Vessel Image')
        axes(hAx(1))
    else
        axes(hAx(2))
    end
    [cROImask{end+1}, cXi{end+1}, cYi{end+1}] = roipoly;
    axes(hAx(1))
    plot(cXi{end}, cYi{end}, 'o-', 'color', mCols(length(cROImask), :))
    axes(hAx(2))
    plot(cXi{end}, cYi{end}, 'o-', 'color', mCols(length(cROImask), :))

    % Check if ESC was pressed
    if isempty(cROImask{end})
        break
    else
        sAns = questdlg('Do want to draw one more ROI?', 'ROI', ...
            'Yes', 'No', 'No');
        if isempty(sAns) || strcmp(sAns, 'No')
            break
        end
    end
end

% Store results
if ~isempty(cROImask{:})
    ISIdata.analysisSignalROI = struct('mROI', cROImask, 'vXi', cXi, 'vYi', cYi);
end

% Close figure
%if ishandle(hFig)
%    close(hFig);
%end

return