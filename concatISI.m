
%pathdir='C:\Documents and Settings\labview_2011\Desktop\20110215_IOS';
% pathdir='C:\Users\Mike\Desktop\dklab\Mike IOS data\2011.02.04';
pathdir='F:\dklab\Mike IOS data\2011.02.16';

cd(pathdir);
d = dir('*.dat');
filenames=struct2cell(d);
filenames=filenames(1,:); %clip everything but the names

% will match first part of filename string, make sure it is exclusive for
%   the file series you want!
filetag='allwhiskers'; 
outfile=[filetag];

filematch=strmatch(filetag,filenames);

%%%%%%%%%%%%%%%
% fixed params
numfullframe=220;
nbinframe=22; %in frames

starttime=uint32(zeros(4,1)); %read as ulong; system dependent, 32 bits on 32bit (i386), 64 bits on x64
size_x=int16(512);
size_y=int16(512);
frame_rate=int16(20);
bin_duration=int16(numfullframe/nbinframe); %should be 10
nsec=uint16(11);
bit_depth1=int16(16); %actually 12 from the camera, but after the conversions it's in 16
ntrials=int32(length(filematch));
intcheck=int32(-9999);

%%%%%%%%%%%%%%%


%all writes should be using little-endian
fout=fopen(outfile,'w');
fwrite(fout,starttime,'uint32');
fwrite(fout,size_x,'int16');
fwrite(fout,size_y,'int16');
fwrite(fout,frame_rate,'int16');
fwrite(fout,bin_duration,'int16');
fwrite(fout,nsec,'uint16');
fwrite(fout,bit_depth1,'int16');
fwrite(fout,ntrials,'int32');
fwrite(fout,zeros(25-4,1),'char'); %skips forward 25, minus 4 bytes for ntrials


fprintf(1,'%s, trial ',filetag);
for fi=1:length(filematch)
    fprintf(1,'%i ',fi-1);
    fid = fopen(filenames{filematch(fi)},'r');
    
    if fid==-1, error('invalid file'); end;

    fseek(fid,12,0);%initial offset
    a = fread(fid,'*int16','b'); %read full file, IN BIG-ENDIAN!!!
    maxp=max(a); minp=min(a);
    frames=reshape(a,512,512,[]);
    fclose(fid);
    
    framestack=cell(nbinframe,1);
    framecount=1;
    for f=1:bin_duration:size(frames,3)
        framestack{framecount}=mean(frames(:,:,f:(f+bin_duration-1)),3);
        framecount=framecount+1;
    end
    
    %trial header
    fwrite(fout,intcheck,'int32');
    fwrite(fout,fi-1,'uint32'); %trial number, starts with zero
    fwrite(fout,-1,'int32'); %stim number, not used
    fwrite(fout,'##:##:##.###','int16'); %trial_times; array of # to fill time chars, 24 bytes
    fwrite(fout,zeros(25,1),'int8'); %skips forward 25 bytes, plus 12, to compensate for incorrect byte size on trial_times
    
    %write the averaged frames
    for f=1:nbinframe
        fwrite(fout,framestack{f},'int16',0,'l');
        if f~=nbinframe %not the last frame
            fwrite(fout,zeros(8,1),'int8'); %skips forward 8; not clear why the file does this between trials
        end
    end
    
    
    
end 
fprintf(1,'\n');

fclose(fout);


