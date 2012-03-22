function ISIdata = ISI_selectROI(ISIdata,prmtsCurrent)

%display reference image and let user select ROI for analysis.
h2RoiSelect = figure('Name','Select ROI');
imRef = imread(fullfile(prmtsCurrent.path2dir,prmtsCurrent.refImage));
[ROImask,xi,yi] = roipoly(imRef);
ISIdata.analysisROI = struct('ROI',ROImask,'xi',xi,'yi',yi);
close(h2RoiSelect);
drawnow

return