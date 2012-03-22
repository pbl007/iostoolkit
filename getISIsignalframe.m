function signalFrame=getISIsignalframe(ISIdata,stimi,figureflag)
%signalFrame=getISIsignalframe(ISIdata,stiminterval,figureflag)
%  given ISIdata, isolate signal frame within given stiminterval.
%  plot the stimulus frame and the signal frame if flag is true

defaultstimi=[0.5 1.5];
if nargin==1
    figureflag=1; %default to plot figure
    stimi=defaultstimi;
elseif nargin==2
    stimi=defaultstimi;
elseif nargin==3
    %do nothing
else
    error('incorrect number of input arguments');
end

% make sure stimi is rounded down to valid frame timepoints
% e.g.  if frame bins are 0.2 sec, then [0.5 1.5] should be
%       rounded to [0.4 1.4]
binpersec = ISIdata.frame_rate/ISIdata.bin_duration;
stimi = floor(stimi ./ (1/binpersec)) .* (1/binpersec); % added, Per 012312

% average deltaR across baseline, stim, and "poststim",

stimFrameix=(ISIdata.nPreStimFrames+1):(ISIdata.nPreStimFrames+ISIdata.nStimFrames);
stimFrame=mean(ISIdata.deltaSignal(:,:,stimFrameix),3);
maxval1=max(max(stimFrame)); minval1=min(min(stimFrame));

% calculate the frames for x1 - x2 sec after stim
%binpersec=ISIdata.frame_rate/ISIdata.bin_duration; % moved up, Per 013212
signalFrameix=(ISIdata.nPreStimFrames+binpersec*stimi(1)):(ISIdata.nPreStimFrames+binpersec*stimi(2));
signalFrame=mean(ISIdata.deltaSignal(:,:,signalFrameix),3);
maxval2=max(max(signalFrame)); minval2=min(min(signalFrame));
clim=[min(minval1, minval2) max(maxval1,maxval2)];

if figureflag
    figure; colormap hot;
    subplot(1,2,1); 
    imagesc(stimFrame); axis image; axis off; set(gca,'clim',clim);
    title('stim frame(s)');
    subplot(1,2,2);
    imagesc(signalFrame); axis image; axis off; set(gca,'clim',clim); 
    title([ num2str(stimi(1)) '-' num2str(stimi(2)) 's poststim']);
end

return