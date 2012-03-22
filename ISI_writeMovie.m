function ISIdata = ISI_writeMovie(ISIdata,prmts)
%Generate movie of averaged frames across trials
% Assumes that ISIdata.deltaSignal exists, with each frame as the averaged
% intrinsic signal across trials. We use the already calculated 
% ISIdata.climAll to scale the image appropriately.
%This is in place of ISI_averageTrialsMovingAverage()
% and ISI_trialsMovingAverage(), which repeat analysis and do so in a
% different manner than what is done in ISI_calc_dRR
%
% Original code was written by Patrick Drew/Pablo Blinder
%modified by Mike Pesavento, 2011.11.11

if ~isfield(ISIdata,'deltaSignal')
    error('deltaSignal does not exist');
end

hmovie=figure('color',[.7 .7 .7],'name',['Trial averaged moving average ' prmts.name]);
set(hmovie,'DoubleBuffer','on');
set(gca,'xlim',[1 ISIdata.frameSizeYX(2)],'ylim',[1 ISIdata.frameSizeYX(1)],...
    'NextPlot','replace','Visible','off');
colormap(gray);
mvi = 0;%frame counter for movie
clear MOVIE ;
secperbin=ISIdata.bin_duration/ISIdata.frame_rate;

%mjp 2010.12.16 fixes crash on filename
filename=fullfile(prmts.path2dir,prmts.name);
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


% % vidobj=videowriter(filename,'fps',10,'width',ISIdata.frameSizeYX(2),'height',ISIdata.frameSizeYX(1));

for fi = 2 : ISIdata.nFramesPerTrial
    img=ISIdata.deltaSignal(:,:,fi);
    imagesc( img);
    
    %title(num2str(fi));
    text(4,ISIdata.frameSizeYX(1)-10, [num2str(fi*secperbin) ' s'],'color','blue','fontweight','bold')
    set(gca,'clim',ISIdata.climAll); colorbar ('horizontal');axis image;axis off
    if fi > ISIdata.nPreStimFrames  && fi <= (ISIdata.nPreStimFrames+ISIdata.nStimFrames)
%         h2rect = rectangle('Position',[10 10 15 15]);
%         set(h2rect,'faceColor','white');
        w=floor(ISIdata.frameSizeYX(1)*0.05); %make patch width 2% of frame
        patch([10 10 10+w 10+w],[10 10+w 10+w 10],'red');
    end
    
    mvi = mvi + 1;
    MOVIE(mvi) = getframe(gca);
%     if mvi>1 && size(MOVIE(mvi).cdata,1)<size(MOVIE(mvi-1).cdata,1)
%         newframe=zeros(size(MOVIE(mvi-1).cdata));
%         newframe(1:size(MOVIE(mvi).cdata,1),1:size(MOVIE(mvi).cdata,2),:)=MOVIE(mvi).cdata;
%         MOVIE(mvi).cdata=newframe;
%     end
    pause(.1)
end

try
%    avimov = avifile(filename,'fps',10,'compression','none','quality',100);
    avimov = avifile(filename,'fps',10,'quality',100);
    avimov = addframe(avimov,MOVIE);
    avimov = close(avimov);
catch ME
    avimov = close(avimov);
    error(ME.message);
end


% movie2avi(MOVIE,filename,'FPS',10);

% movie2avi(MOVIE,filename,'FPS',10,'compression','none');
% movie2avi(MOVIE,filename,'FPS',10,'compression','Cinepak');
