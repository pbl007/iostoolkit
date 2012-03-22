function [contourM ROI]=getbarrelROI(signalFrame,contourquantile)
%assume image has already been masked and filtered.


%first find original quantiles
contourvals=quantile(signalFrame(:),contourquantile); %find histogram quantiles
frameSizeYX(1)=size(signalFrame,2); %Y
frameSizeYX(2)=size(signalFrame,1); %X


himg=imagesc(signalFrame);     axis image; hold on; colormap hot; colorbar
[contourM h]=contour(1:frameSizeYX(2), 1:frameSizeYX(1),...
    signalFrame, contourvals,'color','b');


title('select desired ROI...');

he=imrect(gca);

mask=createMask(he,himg);
ROI=round(getPosition(he));
delete(he);


%maskimg=mask.*croppedimg;
%imagesc(maskimg); axis image; colormap hot; colorbar
% [contourM h]=contour(1:(frameSizeYX(2)-2*cropoffset), 1:(frameSizeYX(1)-2*cropoffset),...
%     maskimg, contourvals,'color','b');
% close(gcf);

if ROI(2)<1, ROI(2)=1; end
if ROI(1)<1, ROI(1)=1; end
if (ROI(2)+ROI(4)) > frameSizeYX(2), ROI(4)=frameSizeYX(2)-ROI(2); end
if (ROI(1)+ROI(3)) > frameSizeYX(1), ROI(3)=frameSizeYX(1)-ROI(1); end


croppedimg=signalFrame(ROI(2):ROI(2)+ROI(4),ROI(1):ROI(1)+ROI(3));
%contourvals=quantile(croppedimg(:),contourquantile); %find histogram quantiles
delete(gca);
himg=imagesc(croppedimg); axis image; colormap hot; colorbar
hold on
[contourM h]=contour(croppedimg, contourvals,'color','b');
close(gcf);

