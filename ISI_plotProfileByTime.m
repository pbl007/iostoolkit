function ISI_plotProfileByTime(ISIdata, prmts, cmap)
% ISI_plotProfileByTime(ISIdata)
% Plot intensities of selected line profiles across frames

% Get/generate ROI
if isfield(ISIdata, 'analysisSignalProfile')
    tProfile = ISIdata.analysisSignalProfile;
else
    tProfile = struct([]);
end

hFig = figure;
hFig2 = figure;

f = size(ISIdata.deltaSignal, 3);
vTime = (((ISIdata.frame_rate/ISIdata.bin_duration):f+((ISIdata.frame_rate/ISIdata.bin_duration)-1)) ...
    / (ISIdata.frame_rate/ISIdata.bin_duration)) ...
    - (ISIdata.nPreStimFrames / (ISIdata.frame_rate/ISIdata.bin_duration)) ...
    - 1;


% Iterate over profiles
for r = 1:length(tProfile)
    figure(hFig)
    hAx = subplot(length(tProfile), 1, r);
    vXi = double(tProfile(r).vXi);
    vYi = double(tProfile(r).vYi);
    
    % Linearly interpolate line segment for each pixel crossed
    vX = min(round(vXi)):max(round(vXi));
    vY = round(interp1(vXi, vYi, vX));
    vXY = unique([vX' vY'], 'rows');

    vY = min(round(vYi)):max(round(vYi));
    vX = round(interp1(vYi, vXi, vY));
    vXY = unique([vXY; vX' vY'], 'rows');

    vXY(any(isnan(vXY), 2), :) = [];
    
    % Initialize smoothing filter
    if ~isnan(prmts.smoothSigma)
        mWin = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
    else
        mWin = NaN;
    end
    
    % Colormap for lines
    mCols = colormap(jet); 
    n = size(ISIdata.deltaSignal, 3); % size of new color map
    m = size(mCols,1);
    t0 = linspace(0,1,m)';
    t = linspace(0,1,n)';
    r = interp1(t0,mCols(:,1),t);
    g = interp1(t0,mCols(:,2),t);
    b = interp1(t0,mCols(:,3),t);
    mCols = [r,g,b];
    colormap(mCols)
    
    % Iterate across frames and compute mean intensity in profile
    mProfileMean = [];
    vWidths = [];
    for f = 1:size(ISIdata.deltaSignal, 3)
        vCol = mCols(f,:);
        
        % Frame signal and error
        mFrameSignal = ISIdata.deltaSignal(:,:,f);
        mFrameError = ISIdata.deltaError(:,:,f);
        if ~isnan(mWin)
            mFrameSignal = single(filter2(mWin, mFrameSignal, 'same'));
            mFrameError = single(filter2(mWin, mFrameError, 'same'));
        end
        
        % Get intensity profile
        vProfileMean = mFrameSignal(sub2ind(size(mFrameSignal), vXY(:,2), vXY(:,1)));
        vProfileError = mFrameError(sub2ind(size(mFrameError), vXY(:,1), vXY(:,2)));
        
        % Get width at half-peak
        [nMin, nMinxIndx] = min(vProfileMean);
        vWidths(f) = max(diff(find(vProfileMean > nMin/2)));
        
        figure(hFig)
        axes(hAx)
        hold on
        plot(vProfileMean, 'color', vCol, 'linewidth', 2)
        mProfileMean = [mProfileMean vProfileMean];
        %errorbar(1:length(vProfileMean), vProfileMean, vProfileError, 'color', mCols(f,:))
        
    end
    
    figure(hFig)
    xlabel('Profile Distance (px)')
    ylabel(sprintf('dF/F (%s)', char(8240)))
    grid on
    axis tight
    hLeg = legend({'Stim' 'Mean' 'Std'}, 'Location', 'Best');
    legend boxoff
    title(sprintf('Intensity Profile - %s : %s', prmts.name, prmts.Whisker{1}), 'interpreter', 'none')
    hAx = colorbar;
    ylabel(hAx, 'Frame')

    
    % Get peaks (location and amplitude)
    [vMin, vMinIndx] = min(mProfileMean, [], 1);
    
    figure(hFig2)
    hAx = subplot(3,1,1);
    hold on
    plot((1:length(vMin))-ISIdata.nPreStimFrames, vMinIndx, '.-');
    xlabel('Post-stim time (frames)')
    ylabel('Peak location (px)')
    axis tight

    % Todo: add quantiles to estimate significance
    hAx = subplot(3,1,2);
    hold on
    plot((1:length(vMin))-ISIdata.nPreStimFrames, vMin, '.-');

    % Plot quantiles (to estimate significance)
    vQuant = quantile(double(ISIdata.deltaSignal(:)), [.5 .75 .95]);
    plot([1 length(vMin)]-ISIdata.nPreStimFrames, [-vQuant(1) -vQuant(1)], 'g')
    text(0, -vQuant(1), '.5')
    plot([1 length(vMin)]-ISIdata.nPreStimFrames, [-vQuant(2) -vQuant(2)], 'g')
    text(0, -vQuant(2), '.75')
    plot([1 length(vMin)]-ISIdata.nPreStimFrames, [-vQuant(3) -vQuant(3)], 'g')
    text(0, -vQuant(3), '.95')
    
    xlabel('Post-stim time (frames)')
    ylabel('Amplitude at peak location (dF/F)')
    set(gca, 'ylim', [min(vMin) max([0 max(vMin)])])
    axis tight

    
    % Plot widths at half-max
    hAx = subplot(3,1,3);
    hold on
    plot((1:length(vMin))-ISIdata.nPreStimFrames, vWidths, '.-');
    xlabel('Frame')
    ylabel('Width at half-peak (px)')
    axis tight
    
end



keyboard

return



function [vC, vInd] = nearest(vA, vB)
% NEAREST Substitute for nearest matrching number
% [C, I] = NEAREST(A,B) substitutes numbers in vector A with the nearest number
% in B and returns the result to C, where the length of C is equal to the
% length of A. The substitution indices of B are returned in I.
%
% Written by Per Magne Knutsen, January 2005

% Reshape vectors
vA = vA(:);
vB = vB(:)';

% Subtract every number in vB from every number in vA
mB = repmat(vB, length(vA), 1);
mA = repmat(vA, 1, length(vB));
mR = mA - mB;
[vMin, vInd] = min(abs(mR),[],2);
vC = vB(vInd);

return
