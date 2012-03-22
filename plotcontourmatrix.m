function c=plotcontourmatrix(C,offset,color,fill,alpha)
% c=plotcontourmatrix(C,offset,color,fill,alpha) plot contour lines on current figure
%
% returns c as a struct array with fields x and y
% C = contour matrix, as returned from using contour() or its lower level derivatives
% offset = pixels to shift contours by if they were obtained from a cropped
%           form of the base image. default=[0 0]
% color = colorspec for line or fill
% fill = flag, if set to 1, use patch to fill the contours; if 0, plot
%           lines. default=0
% alpha = transparency value for patch; 0=transparent, 1=opaque

if nargin ==1
    color='k';
    offset=[0 0];
    fill=0;
    alpha=1;
elseif nargin==2
    color='k';
    fill=0;
    alpha=1;
elseif nargin==3
    fill=0;
    alpha=1;
elseif nargin==4
    alpha=1;
elseif nargin==5
    %should be fine
else
    error('incorrect number of input arguments');
end

%%
c(1).x=[];
c(1).y=[];
curix=1;
count=1;
while curix<=size(C,2)
    nC=C(2,curix);
    curix=curix+1; %shift to start of array
    c(count).x=C(1,curix:(nC+curix-1))+offset(1);
    c(count).y=C(2,curix:(nC+curix-1))+offset(2);
    count=count+1;
    curix=curix+nC;
end


%%
if fill 
    for k=1:length(c)
        patch(c(k).x, c(k).y,color, 'edgealpha',alpha,'facealpha',alpha);
    end
else %fill==0
    for k=1:length(c)
        line(c(k).x, c(k).y,'color',color);
    end
end
