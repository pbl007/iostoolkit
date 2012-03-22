function imgmask=maskvessels(filename, thresh)
% imgmask=maskvessels(filename, thresh)
% create mask of vasculature from image filename, 
%  using default threshold of 0.5 if none passed in

if nargin==1
    thresh=0.5;
elseif nargin==2
    %do nothing
else
    error('incorrect number of arguments');
end
    

%load file
img=imread(filename);
Igray=mat2gray(img); %don't really need to normalize, but do it for clarity
imgmask=im2bw(Igray,thresh);

