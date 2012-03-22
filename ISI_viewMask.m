function ISI_viewMask(filename,thresh)
% view the mask on top of the original vessel file
% used to compare actual coverage, see what areas are being masked out

if nargin==0
    vesselfile=fullfile(basepath,subpath,filename);
    thresh=0.45;
else
    vesselfile=filename;
end
 
vthresh=thresh;

vessel=imread(vesselfile);
vessel=mat2gray(vessel); %don't really need to normalize, but do it for clarity
vmask=im2bw(vessel,vthresh);
vmaskrgb=repmat(vmask,[1 1 3]);
vmaskrgb(:,:,2:3)=0;

figure
hv = imshow(vessel);  
axis image; hold on;
hmask=imshow(vmask);
set(hmask,'alphadata',0.3,'CData',vmaskrgb);

title(['mask threshold = ' num2str(thresh)]);
