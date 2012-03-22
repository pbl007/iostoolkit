function ISIdata = ISI_calc_dRR(ISIdata,prmts,Rnolight)
% ISI_calc_dRR Compute mean and standard deviation of evoked IOS signal
%
%The goal here is to take the difference between the given frame and the
%average baseline, and normalize by the difference between the baseline and
%the nolight condition
%
% (R_i-R_base)/(R_base-R_nl) = signal at frame i
% 
% also, can drop the no light (R_nl) condition if we don't have the file,
% giving
% (R_i-R_base)/(R_base) == R_i/R_base - 1
% This last formula matches that used in 
%      (Kohn, Metz, Quibera, Tommerdahl, and Whitsel, 2000, Neuroscience)
% The main difference with other publications is whether we are blocking
% all stimulus frames (0.5 to 2s after stimulation) together, or separating
% them out across all frames

% Based on ISI_averageTrials() by Patrick Drew and Pablo Blinder
%   Modified by Mike Pesavento, Feb 19 2011
%   Modified by Per M Knutsen, Feb 1 2012
%

if isempty(ISIdata.frameStack)
    error('ISIdata.frameStack is empty or in an old format. Check your input file.');
end

% average trials
trials2use = prmts.Trials2Use;
if isempty(trials2use);
    trials2use = 1:ISIdata.ntrials;
end
useNtrials = length(trials2use);

%initialize sum arrays
deltaSignalSum = cell(1,useNtrials);
deltaErrorSum = cell(1,useNtrials);
for i=1:ISIdata.nFramesPerTrial
    deltaSignalSum{i} = single(zeros(ISIdata.frameSizeYX));
    deltaErrorSum{i} = single(zeros(ISIdata.frameSizeYX));
end

hWait = waitbar(0,'Averaging frames across trials...');
centerfig(hWait, findobj('Tag', 'ISIanalysisGUI_fig'));
drawnow;

% Iterate over trials
normalizingFrames = {};
for ti = trials2use
    if isempty(find(ti==trials2use, 1)) % skip trials we're not using
        continue;
    end
    
    % Get average of baseline for trial ti
    baseframeix = 1:ISIdata.nPreStimFrames;
    baselineFrames = single(cell2mat(ISIdata.frameStack(ti,baseframeix)));
    baselineFrames = reshape(baselineFrames(:),[ISIdata.frameSizeYX ISIdata.nPreStimFrames ]);
    baselineTrial = mean(baselineFrames,3)'; 
    baselineTrialMedianCorr = baselineTrial - median(baselineTrial(:));
    ISIdata.baselineTrial{ti} = single(baselineTrial);
    if ~isempty(Rnolight)
        normalizingFrame = baselineTrial - Rnolight;
    else
        normalizingFrame = baselineTrial;
    end
    normalizingFrames{ti} = normalizingFrame;

    % Add trial signal to deltaSignalSum
    for fi = 1:ISIdata.nFramesPerTrial
        if isempty(ISIdata.frameStack{ti,fi})
            warning('Missing frames in ISIdata.frameStack')
        end
        deltaSignalSum{fi} = deltaSignalSum{fi} + single((ISIdata.frameStack{ti,fi}' - baselineTrial) ./ normalizingFrame);
    end

    % Don't think we actually need to be saving these
    %ISIdata.baselineTrial{ti}=baselineTrial;
    %ISIdata.baselineTrialMedianCorr{ti}=baselineTrialMedianCorr;
    waitbar(ti/max(trials2use), hWait)
end

% Calculate average signal
deltaSignal = single(nan([ISIdata.frameSizeYX ISIdata.nFramesPerTrial]));
deltaSignalMedianCorrected = deltaSignal;
for fi = 1:ISIdata.nFramesPerTrial
    deltaSignal(:,:,fi) = deltaSignalSum{fi}./numel(trials2use); % take manual average
    deltaSignalMedianCorrected(:,:,fi) = deltaSignal(:,:,fi) - median(reshape(deltaSignal(:,:,fi),1,[]));
end
maxp  = max(deltaSignal(:));
minp  = min(deltaSignal(:));
maxpm = max(deltaSignalMedianCorrected(:));
minpm = min(deltaSignalMedianCorrected(:));

% Calculate standard deviation of signal
waitbar(ti/max(trials2use), hWait, 'Computing errors across trials...')
for ti = trials2use
    if isempty(find(ti == trials2use, 1)) % skip trials we're not using
        continue;
    end
    normalizingFrame = normalizingFrames{ti};

    % Add trial signal to deltaSignalSum
    for fi = 1:ISIdata.nFramesPerTrial
        if isempty(ISIdata.frameStack{ti,fi})
            warning('Missing frames in ISIdata.frameStack')
        end
        deltaErrorSum{fi} = (single((ISIdata.frameStack{ti,fi}' - baselineTrial) ./ normalizingFrame) - deltaSignal(:,:,fi)) .^ 2;
    end

    waitbar(ti/max(trials2use), hWait)
end
close(hWait)

% Assign variables to output structure
ISIdata.normFrame = normalizingFrame;
ISIdata.deltaSignal = deltaSignalMedianCorrected;
ISIdata.climAll = [minpm maxpm]; % use values from median corrected r

% Clear some variables to conserve memory
clear baselineFrames baselineTrial baselineTrialMedianCorr deltaSignalSum normalizingFrames normalizingFrame deltaSignalMedianCorrected

% This part may fail with Not Enough Memory. If so, just skip it.
try
    ISIdata.deltaError = single(nan([ISIdata.frameSizeYX ISIdata.nFramesPerTrial]));
    for fi = 1:ISIdata.nFramesPerTrial
        ISIdata.deltaError(:,:,fi) = single( sqrt(deltaErrorSum{fi} ./ numel(trials2use)) ); % standard deviation
    end
catch
    ISIdata.deltaError = [];
end

return
