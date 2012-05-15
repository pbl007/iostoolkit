function ISIdata = ISI_writeMovie(ISIdata,prmts)
% Generate movie of averaged frames across trials
% Assumes that ISIdata.deltaSignal exists, with each frame as the averaged
% intrinsic signal across trials. We use the already calculated 
% ISIdata.climAll to scale the image appropriately.
%
% This is in place of ISI_averageTrialsMovingAverage()
% and ISI_trialsMovingAverage(), which repeat analysis and do so in a
% different manner than what is done in ISI_calc_dRR
%

% Changelog:
%  11/11/11 Adapted from original code by Patrick Drew & Pablo Blinder. MJP
%  03/23/12 Fixed saving AVIs, which was broken. PMK
%  03/23/12 Added option to smooth frames. PMK
%  03/23/12 Frames are now saved in true size. PMK

if ~isfield(ISIdata,'deltaSignal')
    error('deltaSignal does not exist');
end

hmovie = figure('color',[.7 .7 .7],'name',['Trial averaged moving average ' prmts.name]);
set(hmovie,'DoubleBuffer','on');
set(gca,'xlim',[1 ISIdata.frameSizeYX(2)],'ylim',[1 ISIdata.frameSizeYX(1)],...
    'NextPlot','replace','Visible','off');
colormap(gray);
mvi = 0; % frame counter for movie
clear MOVIE ;
secperbin=ISIdata.bin_duration/ISIdata.frame_rate;

% mjp 2010.12.16 fixes crash on filename
filename = fullfile(prmts.path2dir,prmts.name);
[path,name,ext] = fileparts(filename);
if ~strcmpi(ext,'.avi')
    k=strfind(filename,ext);
    if ~isempty(k)
        filename=strcat(filename(1:k(end)-1),'.avi');
    end
end
if isempty(ext)
    filename = strcat(filename,'.avi');
end

mSize = [];
for fi = 1:ISIdata.nFramesPerTrial
    img = ISIdata.deltaSignal(:,:,fi);
    
    % Smooth image
    if ~isnan(prmts.smoothSigma)
        mWin = fspecial('gaussian', prmts.smoothSigma*3, prmts.smoothSigma);
        img = single(filter2(mWin, img, 'same'));
    end
    
    % Mask frame
    if prmts.useManualMask
        tMask = prmts.manualMask;
        if ~isempty(tMask)
            img(~tMask.mROI) = NaN;
            % Crop image
            [vI, vJ] = find(~isnan(img));
            vI = [min(vI) max(vI)];
            vJ = [min(vJ) max(vJ)];
            img = imcrop(img, [vJ(1) vI(1) diff(vJ) diff(vI)]);
        end
    end
    
    % Replace NaN's with zeros
    img(isnan(img)) = 0;
    
    if ~ishandle(hmovie), return, end
    figure(hmovie)
    imagesc(img);
    
    title(sprintf('Frame # %d', fi));
    text(4,ISIdata.frameSizeYX(1)-10, sprintf('%5.2f s', fi*secperbin),'color','k','backgroundcolor','w','fontweight','bold')
    set(gca,'clim',ISIdata.climAll); colorbar ('horizontal');
    if fi > ISIdata.nPreStimFrames  && fi <= (ISIdata.nPreStimFrames+ISIdata.nStimFrames)
        w = floor(ISIdata.frameSizeYX(1)*0.02); % make patch width 2% of frame
        patch([10 10 10+w 10+w],[10 10+w 10+w 10], 'red');
    end
    axis image; axis off
    truesize(gcf, size(img))
    
    mvi = mvi + 1;
    MOVIE(mvi) = getframe(gca);
    mSize(end+1, :) = size(MOVIE(mvi).cdata);
    drawnow
end

% Crop all frames in MOVIE to smallest frame (may vary due to buggy matlab plotting)
vMinSize = min(mSize, [], 1);
for f  = 1:length(MOVIE)
    MOVIE(f).cdata = MOVIE(f).cdata(1:vMinSize(1), 1:vMinSize(2), 1:1:vMinSize(3));
end

try
    movie2avi(MOVIE, filename)
catch ME
    avimov = close(avimov);
    error(ME.message);
end

return