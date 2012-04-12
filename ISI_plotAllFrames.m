function h = ISI_plotAllFrames(ISIdata, clim, prmts, cmap)
% h = ISI_plotAllFrames(ISIdata, clim, colormap)
% Plot grid array of frames over time with timestamps.
%

% Modified by Per M Knutsen , Feb 02 2012
%

h = figure;
set(h, 'color', 'k')
colormap(cmap);

% 1st derivative (change in pixel values over time)
nDer = 0;
if nDer > 0
    mFrames = diff(ISIdata.deltaSignal, nDer, 3);
else
    mFrames = ISIdata.deltaSignal;
end

totalNframes = ISIdata.nFramesPerTrial - nDer;
maxrow = 5; % more than this really won't fit on screen
nrow = maxrow;
ncol = double(ceil(totalNframes/maxrow));
if ncol <= maxrow - 2 % narrow number of cols, try recalculating
    nrow = maxrow - 1;
    ncol = ceil(totalNframes/nrow);
end
nframes = 1;

if ~isnan(prmts.smoothSigma)
    mWin = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
else
    mWin = NaN;
end

for j = 1:nrow
    for i = 1:ncol
        if nframes > totalNframes
            break;
        end
        f = (j-1)*ncol+i;
        
        % Generate axes
        axes('position', [(i*(1/ncol))-(1/ncol) 1-(j*(1/nrow))  1/ncol 1/nrow]) % modified Per Jan 10th 2012
        frame = mFrames(:,:,f); %ISIdata.deltaSignal(:,:,f);
        if ~isempty(frame)
            % median filter
            if isfield(prmts, 'medFilter')
                frame = medfilt2(frame, prmts.medFilter);
            end
            
            % Smooth frame
            if ~isnan(mWin)
                frame = single(filter2(mWin, frame, 'same'));
            end
            
            % Mask frame
            if prmts.useManualMask
                tMask = prmts.manualMask;
                if ~isempty(tMask)
                    frame(~tMask.mROI) = NaN;
                end
            end
            
            % Display frame
            imagesc(frame);
            axis image; axis off;
            if ~isempty(clim)
                set(gca,'clim', clim);
            end
        end
        
        % Set timestamp
        nT = (f+(ISIdata.frame_rate/ISIdata.bin_duration)-1+nDer) ...
            / (ISIdata.frame_rate/ISIdata.bin_duration) ...
            - (ISIdata.nPreStimFrames/(ISIdata.frame_rate/ISIdata.bin_duration)) ...
            - 1;
        
        %vTime = (((ISIdata.frame_rate/ISIdata.bin_duration):f+((ISIdata.frame_rate/ISIdata.bin_duration)-1)) ...
        %    / (ISIdata.frame_rate/ISIdata.bin_duration)) ...
        %    - (ISIdata.nPreStimFrames / (ISIdata.frame_rate/ISIdata.bin_duration)) ...
        %    - 1;
        
        hTxt = text(15, 30, [num2str(nT) ' s'], ...
            'color', 'w', 'backgroundcolor', 'k'); % modified Per Jan 10th 2012
        set(hTxt, 'fontsize', 8)
        nframes=nframes+1;
        
    end
end

set(h, 'position',[20 50 1200 670]);

return

