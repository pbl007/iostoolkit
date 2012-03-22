function ISIdata = ISI_TrialMovingAverage(ISIdata,prmts)
%Generate sequence of moving average for all trials
%% average trials
trials2use = prmts.Trials2Use;
if isempty(trials2use);trials2use = 1 : ISIdata.ntrials;end
nTrials2Use = numel(trials2use);
ISIdata.frameStackMovingAverage = cell(nTrials2Use,ISIdata.nFramesPerTrial-1);
minDelta = inf;
maxDelta = -inf;

%cycle trials
for tri = 1 : nTrials2Use
    %get baseline from averaging the appropriate number of preStim frames
    %for current trial
    baseLineFrames = cell2mat(ISIdata.frameStack(trials2use(tri),1:ISIdata.nPreStimFrames));
    baseLineFrames = reshape(baseLineFrames(:),[ISIdata.frameSizeYX ISIdata.nPreStimFrames ]);
    baseLine = mean(baseLineFrames,3)';
%     baseLineMedianCorrected = baseLine - median(baseLine(:));
    %% generate moving average



    for fi = 2 : ISIdata.nFramesPerTrial
        windowFrames = cell2mat(ISIdata.frameStack(trials2use(tri),fi-1:fi));
        windowFrames = reshape(windowFrames(:),[ISIdata.frameSizeYX 2 ]);
        ISIdata.frameStackMovingAverage{tri,fi - 1}= mean(windowFrames,3)'./baseLine;
        minDelta = min([minDelta min(ISIdata.frameStackMovingAverage{tri,fi-1}(:))]);
        maxDelta = max([maxDelta max(ISIdata.frameStackMovingAverage{tri,fi-1}(:))]);
%         imagesc(ISIdata.frameStackMovingAverage{tri,fi - 1} -...
%             median(ISIdata.frameStackMovingAverage{tri,fi - 1}(:))+1);
%         title(num2str(fi));pause(.1);
    end

end % computing moving average

%%
figure('color',[.7 .7 .7],'name',['Trial averaged moving average ' prmts.name])
mvi = 0;%frame counter for movie
clear MOVIE ;
for tri = 1 : nTrials2Use
    for fi = 3 : ISIdata.nFramesPerTrial - 1
        %imagesc( ISIdata.frameStackMovingAverage{tri,fi} - median(ISIdata.frameStackMovingAverage{tri,fi}(:)));
        imagesc( ISIdata.frameStackMovingAverage{tri,fi}); %don't do the median correction (mjp 2011.11.11)
        text(4,250,[num2str(tri) ' ' num2str(fi)],'color','y')
        % mean( ISIdata.frameStackMovingAverage{tri,fi}(:) - median(ISIdata.frameStackMovingAverage{tri,fi}(:)))
        set(gca,'clim',[minDelta maxDelta]);colorbar ('horizontal');axis image;axis off;colormap bone
        if fi > ISIdata.nPreStimFrames  && fi <= ISIdata.nStimFrames
            h2rect = rectangle('Position',[10 10 15 15]);
            set(h2rect,'faceColor','red');
        
        end
        mvi = mvi + 1;
        MOVIE(mvi) = getframe;
        pause(.01)
    end
end

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

%finally, write the file.
movie2avi(MOVIE,filename,'FPS',10);

