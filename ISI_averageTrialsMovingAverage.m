function ISIdata = ISI_averageTrialsMovingAverage(ISIdata,prmts)
%Generate moving average of averaged frames across trials
% written by Patrick Drew/Pablo Blinder
%modified by Mike Pesavento, 2011.01.27

trials2use = prmts.Trials2Use;
if isempty(trials2use); trials2use = 1 : ISIdata.ntrials;end

%get baseline from averaging the appropriate number of preStim frames
baseLineFrames = cell2mat(ISIdata.trialAveragedFrames(1:ISIdata.nPreStimFrames));
baseLineFrames = reshape(baseLineFrames(:),[ISIdata.frameSizeYX ISIdata.nPreStimFrames ]);
baseLine = mean(baseLineFrames,3);
baseLineMedianCorrected = baseLine - median(baseLine(:));

%% generate moving average
minDelta = inf;
maxDelta = -inf;

for fi = 2 : ISIdata.nFramesPerTrial
    windowFrames = cell2mat(ISIdata.trialAveragedFrames(fi-1:fi));
    windowFrames = reshape(windowFrames(:),[ISIdata.frameSizeYX 2 ]);
    ISIdata.trialAveragedFramesMovingAverage{fi}= mean(windowFrames,3)./baseLine;
    ISIdata.trialAveragedFramesMovingAverage{fi}= ISIdata.trialAveragedFramesMovingAverage{fi} - ...
        median(ISIdata.trialAveragedFramesMovingAverage{fi}(:));
    minDelta = min([minDelta min(ISIdata.trialAveragedFramesMovingAverage{fi}(:))]);
    maxDelta = max([maxDelta max(ISIdata.trialAveragedFramesMovingAverage{fi}(:))]);
end

%%

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
    img=ISIdata.trialAveragedFramesMovingAverage{fi};
    imagesc( img);
    
    %title(num2str(fi));
    text(4,ISIdata.frameSizeYX(1)-10, [num2str(fi*secperbin) ' s'],'color','blue','fontweight','bold')
    set(gca,'clim',[minDelta maxDelta]); colorbar ('horizontal');axis image;axis off
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
%     avimov = avifile(filename,'fps',10,'compression','none','quality',100);
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
