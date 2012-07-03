function ISI_plotSignalByTime(ISIdata, prmts, cmap)
% ISI_plotSignalByTime(ISIdata)

% Get/generate ROI
if isfield(ISIdata, 'analysisSignalROI')
    tROI = ISIdata.analysisSignalROI;
else
    tROI = struct([]);
end

hFig = figure;

mCols = get(hFig, 'DefaultAxesColorOrder');
mCols = [mCols; mCols ./ 2];

f = size(ISIdata.deltaSignal, 3);
vTime = (((ISIdata.frame_rate/ISIdata.bin_duration):f+((ISIdata.frame_rate/ISIdata.bin_duration)-1)) ...
    / (ISIdata.frame_rate/ISIdata.bin_duration)) ...
    - (ISIdata.nPreStimFrames / (ISIdata.frame_rate/ISIdata.bin_duration)) ...
    - 1;

hSignal = subplot(2,4,[1:3 5:7]);
plot([0 prmts.stimDurSec], [0 0], 'r', 'linewidth', 10)
hold on

hAutoCorr = subplot(2,4, 4); hold on

hWait = waitbar(0,'Computing ROI amplitude across frames...');
centerfig(hWait, hFig);

% Iterate over ROIs
for r = 1:length(tROI)
    mROI = double(tROI(r).mROI);
    
    % Iterate across frames and compute mean intensity in ROI
    vMean = zeros(1, size(ISIdata.deltaSignal, 3));
    vErr = zeros(1, size(ISIdata.deltaSignal, 3));
    mROI(mROI == 0) = NaN;
    
    drawnow;
    
    if ~isnan(prmts.smoothSigma)
        mWin = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
    else
        mWin = NaN;
    end
    
    for f = 1:size(ISIdata.deltaSignal, 3)
        
        % Frame signal
        mFrameSignal = ISIdata.deltaSignal(:,:,f);
        if ~isnan(mWin)
            mFrameSignal = single(filter2(mWin, mFrameSignal, 'same'));
        end
        mFrameSignal = mFrameSignal .* mROI;
        vMean(f) = nanmean(mFrameSignal(:));
        
        % Frame error
        if ~isempty(ISIdata.deltaError)
            mFrameError = ISIdata.deltaError(:,:,f);
            if ~isnan(mWin)
                mFrameError = single(filter2(mWin, mFrameError, 'same'));
            end
            mFrameError = mFrameError .* mROI;
            vErr(f) = nanmean(mFrameError(:)) / sqrt(size(ISIdata.deltaSignal, 3));
        end
        
        waitbar( (f+(size(ISIdata.deltaSignal, 3)*(r-1))) / (size(ISIdata.deltaSignal, 3) * length(tROI)), hWait);
    end
    
    figure(hFig)
    axes(hSignal)
    if ~isempty(ISIdata.deltaError)
        %mean_error_plot(vMean.*1000, vErr.*1000, mCols(r,:), vTime)
        errorbar(vTime, vMean.*1000, vErr.*1000, vErr.*1000, 'color', mCols(r,:))
    end
    hold on
    plot(vTime, vMean.*1000, 'color', 'k')
    xlabel('Time (s)')
    ylabel(sprintf('dF/F (%s)', char(8240)))
    set(gca, 'ylim', ISIdata.climAll.*1000)
    grid on
    axis tight
    hLeg = legend({'Stim' 'Mean' 'Std'}, 'Location', 'Best');
    legend boxoff
    title(sprintf('ROI Signal - %s : %s', prmts.name, prmts.Whisker{1}), 'interpreter', 'none')

    % Auto-correlations
    axes(hAutoCorr)
    vAutoCorr = xcorr(diff(vMean), 'coeff');
    vACorrTime = linspace(-length(vAutoCorr)/2, length(vAutoCorr)/2, length(vAutoCorr));
    nFramesSec = ISIdata.frame_rate/ISIdata.bin_duration;
    vACorrTime = vACorrTime ./ nFramesSec;
    plot(vACorrTime, vAutoCorr, 'color', mCols(r,:))
    set(gca, 'xlim', [-20 20]) % seconds
    xlabel('Time (s)')
    ylabel('C')
    title('Auto-correlation')
    
    % TODO
    % Plot power-spectrum and coherence
end
close(hWait)


return


function [hPlot,hFill] = mean_error_plot(vMean, vError, vColor, vX)
vMean = reshape(vMean, length(vMean), 1);
vError = reshape(vError, length(vError), 1);
if ~exist('vX')
    vXt = (1:length(vError))';
else
    vXt = vX';
end
vXb = flipud(vXt);
vYt = vMean + vError;
vYb = flipud(vMean - vError);
hFill = fill([vXt;vXb], [vYt;vYb], vColor, 'EdgeColor', vColor); hold on;
hPlot = plot(vXt, vMean, 'k', 'LineWidth', 1); hold off
%set(hFill,'FaceAlpha',.5,'EdgeAlpha',.5) % transparency
return;
