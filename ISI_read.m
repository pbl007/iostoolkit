function [ISIdata] = ISI_read(prmts)
% ISI_READ  Load ISI data
% Load trials from binary file
%
% If function reports an Out of Memory error, enable pixel binning to
% reduce memory usage when the function runs.
%

% Mike Pesavento, 2011.03
%  Modified from Pablo's code, engineered to match the output of ISI_recordframes_piezo.vi in LabVIEW
%  all data formats should be in little-endian, but linux may not read the 8 byte spacers correctly
%  also converts int16 data to single, instead of double

% Per M Knutsen, Feb 02 2012
%  Modified to allow pixel binning during file read.

ISIdata = [];
path2file = fullfile(prmts.path2dir,prmts.name);
fid=fopen(path2file);
if fid==-1
    error(['failed to open ' path2file]);
    return;
end

thestarttime = fread(fid,4,'ulong');
size_x = fread(fid, 1, 'int16');
size_y = fread(fid, 1, 'int16');
frame_rate = fread(fid,1,'int16');
bin_duration = fread(fid,1,'int16');
nsec = fread(fid,1,'uint16');
bit_depth1 = fread(fid,1,'int16');
fpos = ftell(fid);
ntrials = fread(fid,1,'int32');
nFramesPerTrial = nsec*(frame_rate/bin_duration);
bin_duration_sec = bin_duration / frame_rate;
if ~isempty(prmts.Trials2Use)
    trial_range = sprintf('[%d %d]', prmts.Trials2Use(1), prmts.Trials2Use(end));
else
    trial_range = sprintf('[0 %d]', ntrials);
end
fprintf('\nFilename:\t\t%s\nTrials:\t\t\t%d\nFrame Rate:\t\t%d frames/s\nBin Duration:\t%d frames / %.2f s\nFrame Size:\t\t%dx%d px\nFrames/trial:\t%d\nBit Depth:\t\t%d\nTrial Duration:\t%d s\nTrials Used:\t%s\n',...
    prmts.name, ntrials, frame_rate, bin_duration, bin_duration_sec, size_x,size_y,nFramesPerTrial,bit_depth1,nsec,trial_range);

if ~prmts.DoLoad
    return
end

fseek(fid, 25-4, 0);

hWait = waitbar(0,'Loading frames. Please wait...');
centerfig(hWait, findobj('Tag', 'ISIanalysisGUI_fig'));
drawnow;

ISIdata.frameStack = cell(ntrials,(nsec*(frame_rate/bin_duration)));
size(ISIdata.frameStack);
for k = 1:ntrials
    %read trial header
    %intcheck = fread(fid,1,'int32','l', sMachineformat);        % should be -9999
    intcheck = fread(fid,1,'int32');        % should be -9999
    if intcheck ~= -9999
        warning('intcheck does not return -9999 in ISI_read.m. Entering debugging mode.')
        keyboard
    end
    trial = fread(fid,1,'uint32');          % trial #, starting with 0
    stimnum = fread(fid,1,'int32');
    trial_times = fread(fid,24,'signed char');   % 24
    fseek(fid, 25, 0);
    fseek(fid, -8, 0); % back up by 4 int16 because not writing array size, like earlier version

    for m = 1:((1/bin_duration_sec)*nsec) % modified Per Jan 10th 2012
        mFrame = single([]);
        mFrame = single(fread(fid, [size_x size_y], prmts.precision));

        % Bin pixels
        if isfield(prmts, 'imgBin')
            if prod(prmts.imgBin) > 1
                mFrame = downsamp2d(mFrame, prmts.imgBin);
            end
        end
        
        % Assign 
        ISIdata.frameStack{k,m} = single([]);
        ISIdata.frameStack{k,m} = mFrame;
    end
    if ishandle(hWait)
        waitbar(k/ntrials, hWait)
    else
        error('File read aborted.');
    end
end
close(hWait)
fclose(fid);

%store parameters relevant for further analysis
ISIdata.ntrials = ntrials;
ISIdata.nsec = nsec;
ISIdata.frame_rate = frame_rate;
ISIdata.bin_duration = bin_duration;
ISIdata.nFramesPerTrial = nFramesPerTrial;

% Frame size
ISIdata.frameSizeYX = [size_y size_x];
if isfield(prmts, 'imgBin')
    if prod(prmts.imgBin) > 1
        ISIdata.frameSizeYX = ISIdata.frameSizeYX ./ prmts.imgBin;
    end
end

ISIdata.nPreStimFrames = (frame_rate./bin_duration) * prmts.preStimDurSec;
ISIdata.nStimFrames = (frame_rate./bin_duration) * prmts.stimDurSec;
ISIdata.nPostStimFrames = nFramesPerTrial - (ISIdata.nPreStimFrames + ISIdata.nStimFrames);

return
