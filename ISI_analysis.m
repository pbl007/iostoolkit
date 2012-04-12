function [prmts, ISIdata] = ISI_analysis(prmts)
% [prmts,ISIdata] = ISI_analysis(prmts)   
% Main function for IOS analysis.


%  modified from Pablo Blinder, ISI_analysis_trialBased.m
%  M.Pesavento 2010.12.30
%
% 2011.07.04 mjp changed from script to function, removed multifile
%           specifications, cleaned up function calls
%
% 2012.02.02 Per Magne Knutsen
%            Change to transition all functions to be called from GUI.
% 2012.04.03 Added manual mask. PMK
%

sColMap = 'gray';

% Set the no light file data, subtracting the DC offset
if prmts.useNoLight
    % Get the nolight data
    if ~exist(fullfile(prmts.Rnolight.path2dir, prmts.Rnolight.name),'file')
        choice=questdlg('Continue without NoLight file?','Continue?','Yes','No','Yes');
        if strcmpi(choice,'No')
            return;
        end
        drawnow
        Rnolight=[];
    else
        darkdata = ISI_read(prmts.Rnolight);
        % make darkfield average
        Rnolightframes=nan([darkdata.frameSizeYX darkdata.nFramesPerTrial]);
        Rnolighttrial=nan([darkdata.frameSizeYX darkdata.ntrials]);
        for i=1:darkdata.ntrials
            for j=1:darkdata.nFramesPerTrial
                Rnolightframes(:,:,j)=cell2mat(darkdata.frameStack(i,j))'; %transpose for correct orientation
            end
            Rnolighttrial(:,:,i) = mean(Rnolightframes,3);
        end
        Rnolight = mean(Rnolighttrial,3);
    end
else
    drawnow
    Rnolight=[];
end

% Cycle through files
for iFile = 1 : numel(prmts.filesQueue)
    % Get data from GUI, if already set
    hGUI = findobj('Tag', 'ISIanalysisGUI_fig');
    tUserData = get(hGUI, 'UserData');
    ISIdata = [];
    if ~isempty(tUserData)
        if isfield(tUserData, 'frameStack')
            ISIdata = tUserData;
        end
    end
    
    % Load data if 1) this is the first time loading this trial, or 2) the
    % Load Data button was pressed (but not the Average Trials or Analyze
    % buttons)
    if isempty(ISIdata) || ...
            (prmts.filesQueue.DoLoad && ~prmts.filesQueue.DoTrialAverage && ~prmts.filesQueue.DoAnalyze)
        % read data
        ISIdata = ISI_read(prmts.filesQueue(iFile));
        
        % store data in GUI
        set(hGUI, 'UserData', ISIdata)
    end

    % Exit here if we should not average trials
    if ~prmts.filesQueue.DoTrialAverage
        return
    end
    
    % Select ROI from vessel map
    %if prmts.filesQueue(iFile).selectROI
    %    ISIdata = ISI_selectROI(ISIdata, prmts.filesQueue(iFile));
    %end

    % Compute trial-averaged frames, if:
    % 1) averages have not already been computed, or
    % 2) the Average Trials button was pressed
    if ~isfield(ISIdata, 'deltaSignal') || ...
            (prmts.filesQueue.DoTrialAverage && ~prmts.filesQueue.DoAnalyze)
        ISIdata = ISI_calc_dRR(ISIdata, prmts.filesQueue(iFile),Rnolight);
        set(hGUI, 'UserData', ISIdata) % store data in GUI
    end

    % Exit now if we should not run analysis
    if ~prmts.filesQueue.DoAnalyze
        return
    end
        
    ISIdata.climAll = prmts.filesQueue(iFile).climAll ./ 1000;
    
    % Get the average frame over stimulus interval, suppressing the comparison figure
    signalFrame = getISIsignalframe(ISIdata, prmts.filesQueue(iFile).stimInterval, 0);
    
    % Plot and save all frames
    if (prmts.saveFigAll)
        ISI_plotAllFrames(ISIdata, ISIdata.climAll, prmts.filesQueue, sColMap);
        drawnow;
        savefilename = fullfile(prmts.filesQueue.path2dir, prmts.filesQueue.name);
        [p n e] = fileparts(savefilename);
        saveas(gcf, fullfile(p,[n '_AllFrames.png']), 'png');
    end
    
    ISIdata.stimInterval = prmts.filesQueue(iFile).stimInterval; %interval we averaged over
    ISIdata.signalFrame = signalFrame;
    
    % Plot trial-average signal of ROI across time
    if prmts.filesQueue(iFile).selectSignalROI
        %ISI_plotContourByTime(ISIdata, prmts.filesQueue(iFile), sColMap);
        ISIdata = ISI_selectSignalROI(ISIdata, prmts.filesQueue(iFile), sColMap);
        if isfield(ISIdata, 'analysisSignalROI')
            ISI_plotSignalByTime(ISIdata, prmts.filesQueue(iFile), sColMap);
        end
    end

    % Plot spatial intensity profile of user-selected line across time
    if prmts.filesQueue(iFile).selectSignalProfile
        ISIdata = ISI_selectSignalProfile(ISIdata, prmts.filesQueue(iFile), sColMap);
        if isfield(ISIdata, 'analysisSignalProfile')
            ISI_plotProfileByTime(ISIdata, prmts.filesQueue(iFile), sColMap);
        end
    end
    
    % Isolate barrel
    if prmts.filesQueue(iFile).runBarrelFinder
        % Create vessel mask
        vesselfile=fullfile(prmts.filesQueue(iFile).path2dir, prmts.filesQueue(iFile).refImage);
        vthresh=prmts.filesQueue(iFile).maskthresh;
        
        if ~prmts.useVesselMask || ~exist(vesselfile,'file')
            vmask = ones(ISIdata.frameSizeYX);
        else
            vmask = imread(vesselfile);
            vmask = mat2gray(vmask); %don't really need to normalize, but do it for clarity
            vmask = im2bw(vmask,vthresh);
        end
        ISIdata.vesselmask = vmask;
        
        % Select contour for barrel
        ISIdata = ISI_isolateBarrel(ISIdata, signalFrame, prmts.filesQueue(iFile),prmts.saveFig);
        drawnow;
    end

    % Generate moving average movie of averaged trials and for all frames
    if prmts.saveMovie
        %create movie of averaged trials
        ISIdata = ISI_writeMovie(ISIdata, prmts.filesQueue(iFile));
        %ISIdata = ISI_averageTrialsMovingAverage(ISIdata, prmts.filesQueue(iFile));
        %creates full movie of all trials
        %ISIdata = ISI_TrialMovingAverage(ISIdata,prmts.filesQueue(iFile));
    end    
    
    
    % save mat file  
    drawnow; % update the figures before saving the mat file
    pause(0.01);
    savefilename = fullfile(prmts.filesQueue(iFile).path2dir,prmts.filesQueue(iFile).name);
    [p n e]=fileparts(savefilename);
    
    if prmts.saveToMat
        if ~strcmp(e,'.mat')
            matfilename = fullfile(p, [n '.mat']);
        end
        
        %cut out the frameStack field, since it's huge and unneccessary
        %mjp 2011.10.13
        ISIdata = rmfield(ISIdata, 'frameStack');
        
        save(matfilename, 'ISIdata');
    end

    % store data and analyzed parameters in GUI
    set(hGUI, 'UserData', ISIdata)

    
end %cycling files