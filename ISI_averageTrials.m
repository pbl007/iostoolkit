function ISIdata = ISI_averageTrials(ISIdata,prmts)

%% average trials
trials2use = prmts.Trials2Use;
if isempty(trials2use);trials2use = 1 : ISIdata.ntrials;end

%raw data average
for fi = 1 : ISIdata.nFramesPerTrial
    %extract
    thisFrame = cell2mat(ISIdata.frameStack(trials2use,fi))'; %image frame gets transposed!
    thisFrame = reshape(thisFrame(:),[ISIdata.frameSizeYX numel(trials2use) ]);
    ISIdata.trialAveragedFrames{fi} = mean(thisFrame,3)';
end

%get baseline from averaging the appropriate number of preStim frames
baseLineFrames = cell2mat(ISIdata.trialAveragedFrames(1:ISIdata.nPreStimFrames));
baseLineFrames = reshape(baseLineFrames(:),[ISIdata.frameSizeYX ISIdata.nPreStimFrames ]);
baseLine = mean(baseLineFrames,3)';
baseLineMedianCorrected = baseLine - median(baseLine(:));

%average stim frame and median correct
stimFrames = cell2mat(ISIdata.trialAveragedFrames((1:ISIdata.nStimFrames) + ISIdata.nPreStimFrames));
stimFrames = reshape(stimFrames(:),[ISIdata.frameSizeYX ISIdata.nStimFrames ]);
stimAvg = mean(stimFrames,3)';
deltaStimAvg = stimAvg./baseLine;
deltaStimAvgMedianCorrected = deltaStimAvg - median(deltaStimAvg(:));

%average post-stim frame and median correct
postStimFrames = cell2mat(ISIdata.trialAveragedFrames((1:ISIdata.nPostStimFrames) + ISIdata.nStimFrames));
postStimFrames = reshape(postStimFrames(:),[ISIdata.frameSizeYX ISIdata.nPostStimFrames ]);
postStimAvg = mean(postStimFrames,3)';
deltaPostStimAvg = postStimAvg./baseLine;
deltaPostStimAvgMedianCorrected = deltaPostStimAvg - median(deltaPostStimAvg(:));

%display and equalize colormap limits
if isfield(prmts,'doDisplayTimeAverage')
    if prmts.doDisplayTimeAverage
        figure('color',[.7 .7 .7],'name',['stim & post-stim avg (median substracted) ' prmts.name],'windowStyle','docked')
        subplot(1,2,1);imagesc(deltaStimAvgMedianCorrected + 1);axis image;axis off;title('Stim');
        subplot(1,2,2);imagesc(deltaPostStimAvgMedianCorrected + 1);axis image;axis off;title('Post-stim');
        h2axes = findobj(gcf,'type','axes');
        clim = cell2mat(get(h2axes,'clim'));
        set(h2axes,'clim',[min(clim(:,1)) max(clim(:,2))])
        colormap copper
        subplot(1,2,1);colorbar horizontal
        subplot(1,2,2);colorbar horizontal
    end
end






%output
ISIdata.baseLine = baseLine;
ISIdata.baseLineMedianCorrected = baseLineMedianCorrected;
ISIdata.stimAvg = stimAvg;
ISIdata.deltaStimAvg = deltaStimAvg;
ISIdata.deltaStimAvgMedianCorrected = deltaStimAvgMedianCorrected;
ISIdata.postStimAvg = postStimAvg;
ISIdata.deltaPostStimAvg = deltaPostStimAvg;
ISIdata.deltaPostStimAvgMedianCorrected = deltaPostStimAvgMedianCorrected;



