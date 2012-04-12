function SVD_analysis(~, ~, handles)


%TODO - move svd dependent parameters up to a user interface

hWait = waitbar(0,'SVD analysis ...');
prmts.hWait = hWait;

%HARDCODED VALUES
prmts.precision='int16';
prmts.DoLoad=1;

prmts.svd.numModesToCompute = []; %leave empty to compute all
prmts.svd.numModesToKeep = 3;
prmts.svd.downsampleFactor = 0.5; % set to 0 to skip, values larger than 1 will be truncated
prmts.svd.DoDenoise = 0;
prmts.svd.DoDisplayEigenvalues = 0;
prmts.svd.DoDisplayMeanSubData = 0;
prmts.svd.DoDisplayDenoised = 0;
prmts.svd.DoDisplaySpatialModes = 0;
prmts.svd.DoSaveSpatialModesToFile = 1;

%validate svd prmts
prmts = SVD_analysis_checkPrmts(prmts);


% Figure out file task(s)
if get(handles.process_all_files, 'value')
    % Get list of all .dat files in directory
    sFileList = dir(fullfile(handles.pathstr,'*.dat'));
else
    sFileList = get(handles.data_filename, 'string');
end

if isempty(sFileList)
    warndlg('No file has been selected.', 'ISI')
    return;
end

% populate prmts structure with values that do not change per each file
prmts.path2dir = handles.pathstr;
prmts.Trials2Use =str2double([get(handles.trials2use_Start,'string') get(handles.trials2use_End,'string')]);
if isnan(prmts.Trials2Use); prmts.Trials2Use = []; end
prmts.preStimDurSec = str2double(get(handles.preStimDurSec,'string'));
prmts.stimDurSec = str2double(get(handles.stimDurSec,'string'));

% set task size
if isstruct(sFileList)
    nLoop = length(sFileList);
else
    nLoop = 1;
    tmp.name = sFileList;
    clear sFileList;
    sFileList = tmp;
end

%create target fig export directories if not existent
eigDirTarget = fullfile(prmts.path2dir,'Eigenvalues');
if ~isdir(eigDirTarget);mkdir(eigDirTarget);end

if prmts.svd.DoSaveSpatialModesToFile
    spatialModesDirTarget = fullfile(prmts.path2dir,'SpatialModes');
    if ~isdir(spatialModesDirTarget); mkdir(spatialModesDirTarget);end
end

% run all task(s)
for iFILE = 1:nLoop
    
    prmts.name = sFileList(iFILE).name;
    [~,prmts.baseName] = fileparts(prmts.name);
    
    waitbar(iFILE/nLoop,prmts.hWait,sprintf('SVD analysis - file %d/%d',iFILE,nLoop))
    
    %get trail averated ISI data (including frameStack) concatenated for svd  (3D->2D, time x space)
    X = SVD_analysis_getData(prmts);
    
    %downsample (if specified) and mean substract
    [X prmts] = SVD_analysis_prepDataForSVD(X,prmts);
    
    %run svd and display/save as required
    waitbar(iFILE/nLoop,prmts.hWait,sprintf('SVD analysis - file %d/%d - computing svd',iFILE,nLoop))
    SVD_analysis_runSVD(X,prmts)
    
    
    
end

close(hWait)



%------------------------------------------- FUNCTIONS ---------------------------------------------------


% prmts = SVD_analysis_checkPrmts(prmts)              validate prmts
% X = SVD_analysis_getData(prmts)
% X = SVD_analysis_prepDataForSVD(X,prmts)
% SVD_analysis_runSVD(X,prmts)

%----------------------------------------------------------------
function X  = SVD_analysis_getData(prmts)
% load data
[ISIdata] = ISI_read(prmts);
% reshape data from 3d to 2d by concatenating images as row-vectors across columns and time across rows
ISIdata = ISI_averageTrials(ISIdata,prmts);
X=cell2mat(ISIdata.trialAveragedFrames);
X=reshape(X(:),[ISIdata.frameSizeYX ISIdata.nFramesPerTrial]); % this reshapes to a 3D time series

%----------------------------------------------------------------
function [X prmts] = SVD_analysis_prepDataForSVD(X,prmts)

[ny,nx,nt] = size(X);

%downsample if required

if prmts.svd.downsampleFactor>0
    for iT = 1 : nt
        tmp(:,:,iT) = imresize(X(:,:,iT),prmts.svd.downsampleFactor);
    end
    X=tmp;clear tmp
    [ny,nx,nt] = size(X);
end
X = reshape(X, [ny * nx  nt])';

%mean substract
meanX = mean(X,1);
X = X - repmat(meanX,[nt 1]);

if prmts.svd.DoDisplayMeanSubData
    %display mean substraced data
    figure('Color','w','Name',sprintf('%s - Reshaped and mean substracted data',prmts.name));
    h2sb1 = subplot(2,1,1);
    imagesc(X);
    colorbar;ylabel('time');xlabel('Image space(x,y) -> s')
    thisXlim = get(gca,'xlim');
    colormap copper
    
    h2sb2 = subplot(2,1,2);
    plot(meanX,'k');set(gca,'xlim',thisXlim);
    hold on;
    windowSize = 50;
    meanXfilt = filter(ones(1,windowSize)/windowSize,1,meanX);
    plot(meanXfilt,'r-','lineWidth',1)
    ylabel('Average Intensity (A.U.)');
    xlabel('Time x trials (A.U.)');
    hc=colorbar;set(hc,'visible','off')
    linkaxes([h2sb1 h2sb2],'x')
end

prmts.nx = nx;
prmts.ny = ny;
prmts.nt = nt;



%----------------------------------------------------------------
function prmts = SVD_analysis_checkPrmts(prmts)
% check analysis settings
numModesToCompute = prmts.svd.numModesToCompute;
downSampleFactor = prmts.svd.downsampleFactor;

if downSampleFactor >1 ; downSampleFactor = 1; end

if prmts.svd.DoDenoise
    numModesToKeep = prmts.svd.numModesToKeep;
    if numModesToKeep > numModesToCompute;
        numModesToCompute = numModesToKeep;
        disp('Increasing number of modes to compute to make it equal to number of modes to keep');
    end
    % if not denoising make sure we do not attemp to plot
    prmts.svd.DoDisplayDenoised = 0;
    
end

prmts.svd.numModesToCompute = numModesToCompute;
prmts.svd.downsampleFactor = downSampleFactor;


%----------------------------------------------------------------
function SVD_analysis_runSVD(X,prmts)

numModesToCompute  = prmts.svd.numModesToCompute;
if isempty(numModesToCompute);numModesToCompute = prmts.nt;end

[u,s,v] = svds(double(X),numModesToCompute);

% display singular values
h2fig = figure('Color','w','Name',sprintf('%s_Eigenvalues',prmts.baseName),'visible','off');
semilogy(diag(s.^2),'ko-.','markerSize',12);
ylabel('Log(\lambda_i^2)')
xlabel('Ranked index (i)')
box off



%save analysis and figure
save(fullfile(prmts.path2dir,[prmts.baseName '_svd_analysis.mat']),'X','u','s','v');
figName = get(h2fig,'name');
export_fig(fullfile(prmts.path2dir,['Eigenvalues' filesep figName]),'-eps','-png');
% saveas(h2fig,fullfile(prmts.path2dir,['Eigenvalues' filesep figName]),'png')
% saveas(h2fig,fullfile(prmts.path2dir,['Eigenvalues' filesep figName]),'eps')
if prmts.svd.DoDisplayEigenvalues;set(h2fig,'visible','on');else close (h2fig);end


% denoise and display
%keep only the first M
if prmts.svd.DoDisplayDenoised
    M = prmts.svd.numModesToKeep;
    Xrecon = u(:,1:M) * s(1:M,1:M) * v(:,1:M)';
    figure('Name',sprintf('%s_Denoised_time_series',prmts.baseName),'color','w');
    subplot(2,1,1)
    title('Original')
    %     imagesc(reshape(X',[nx ny * nt]))
    imagesc(X)
    subplot(2,1,2)
    title(sprintf('Reconstructed with %d',M));
    imagesc(Xrecon)
    colormap copper
end


%% display modes (must create figure to save it
if prmts.svd.DoDisplaySpatialModes || prmts.svd.DoSaveSpatialModesToFile
    
    h2fig = figure('Color','w','Name',sprintf('%s_SpatialModes',prmts.baseName),'visible','on');
    nspcols = ceil(sqrt(numModesToCompute));
    nsprows = ceil(numModesToCompute/nspcols);
    
    eigval = diag(s);
    eigvalProp = eigval./sum(eigval)*100;
    
    for iSP = 1 : numModesToCompute
        subplot(nsprows,nspcols,iSP)
        imagesc(reshape(v(:,iSP),prmts.ny,prmts.nx));axis image;axis off
        title(sprintf('%3.2f',eigvalProp(iSP)))
    end
    
    colormap copper
    
    if prmts.svd.DoSaveSpatialModesToFile
        figName = get(h2fig,'name');
        set(h2fig,'units','norm','pos',[0 0 1 1]);drawnow
        export_fig(fullfile(prmts.path2dir,['SpatialModes' filesep figName]),'-native','-eps','-png');
        %         saveas(h2fig,fullfile(prmts.path2dir,['SpatialModes' filesep figName]),'png')
        %         saveas(h2fig,fullfile(prmts.path2dir,['SpatialModes' filesep figName]),'eps')
    end
    
    if ~prmts.svd.DoDisplaySpatialModes;close(h2fig);else set(h2fig,'visible','on');end
    
    
end

