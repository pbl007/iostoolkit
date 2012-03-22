% accumdiffISI
% takes bins with accumulated data, finds the difference for each trial,
% and dumps to correctly formatted binary file
%
%use this with data from mkII_free_run.vi, which dumps the accumulators to
%a file once for each trial. need to find difference between successive
%accumulators to get raw data

%pathdir='F:\dklab\Mike IOS data\2011.02.16';
% pathdir='E:\IOS data';
pathdir='C:\Users\Mike\Desktop\dklab\Mike IOS data\2011.02.17';

infile='B1_3_213254.dat';
fullinfile=fullfile(pathdir,infile);

outfile=fullfile(pathdir, [infile(1:end-4) '_hdr.dat']); 


%%%%%%%%%%%%%%%
% fixed params
numfullframe=220;
nbinframe=22; %in frames

starttime=uint32(zeros(4,1)); %read as ulong; system dependent, 32 bits on 32bit (i386), 64 bits on x64
size_x=512;
size_y=512;
frame_rate=20;
bin_duration=numfullframe/nbinframe; %should be 10
nsec=11;
bit_depth1=16; %actually 12 from the camera, but after the conversions it's in 16
%ntrials=int32(length(filematch));
ntrials=nan;
intcheck=-9999;

%%%%%%%%%%%%%%%



%all writes should be using little-endian
fout=fopen(outfile,'w+');
if fout==-1
    error('couldn''t open output file');
end
fwrite(fout,starttime,'uint32');
fwrite(fout,size_x,'int16');
fwrite(fout,size_y,'int16');
fwrite(fout,frame_rate,'int16');
fwrite(fout,bin_duration,'int16');
fwrite(fout,nsec,'uint16');
fwrite(fout,bit_depth1,'int16');
ntrialposition=ftell(fout);
fwrite(fout,9999,'int32'); %placeholder for ntrial
%fwrite(fout,ntrials,'int32');
fwrite(fout,zeros(25-4,1),'char'); %skips forward 25, minus 4 bytes for ntrials


fprintf(1,'%s, trial ',infile);

fid = fopen(fullinfile,'r');
if fid==-1, error('invalid file'); end;
%fseek(fid,12,0);%initial offset

%read in all frames from all trials
%should be 22 bins in x trials
nbin=22;
sizex=512;
sizey=512;
np=sizex*sizey;

lasttrial=int32(zeros(sizex*sizey*nbin,1));
trialcount=0;

while ~feof(fid)
    fprintf(1,'%i ',trialcount); %start with trial zero

    trialcount=trialcount+1; %index by 1, not zero
    %read one trial in at a time
    curtrial = fread(fid,sizex*sizey*nbin,'*int32','l');
    if size(curtrial,1)~=size(lasttrial)
        warning('frame data not the same size');
    end
    %write trial header
    fwrite(fout,intcheck,'int32');
    fwrite(fout,trialcount-1,'uint32'); %trial number, starts with zero
    fwrite(fout,-1,'int32'); %stim number, not used
    fwrite(fout,'##:##:##.###','int16'); %trial_times; array of # to fill time chars, 24 bytes
    fwrite(fout,zeros(25,1),'int8'); %skips forward 25 bytes, plus 12, to compensate for incorrect byte size on trial_times

    accumdiff=(curtrial-lasttrial)./bin_duration;
    for fi=1:nbinframe
        %we write out the difference between the trials
        fwrite(fout,accumdiff(np*(fi-1)+(1:np)) ,'int16',0,'l');
        if fi~=nbinframe %not the last frame
            fwrite(fout,zeros(8,1),'int8'); %skips forward 8; not clear why the file does this between trials
        end
    end
    %read next word to see if at eof
    test=fread(fid,1,'int16');
    if ~isempty(test)
        fseek(fid,-2,0); %go back two bytes, or one int16
    elseif feof(fid)        
        break; %get out!
    else
        disp('dunno what happened');
    end
    lasttrial=curtrial;
end
fclose(fid);

fseek(fout,ntrialposition,'bof');
if fread(fout,1,'int32')==9999
    fseek(fout,-4,0);%backup one int32
    fwrite(fout,trialcount,'int32'); %and write
else
    warning('incorrect file position for trialcount');
end

fprintf(1,'\n');

stat=fclose(fout);
if stat==-1, error('output file didn''t close properly'); end



