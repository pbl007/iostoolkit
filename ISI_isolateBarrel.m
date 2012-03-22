function ISIdata=ISI_isolateBarrel(ISIdata,signalFrame,prmtsCurrent,saveFig)
% ISIdata=isolatebarrel(ISIdata,signalFrame,contourquantile,selectROI)
% given signalFrame and desired quantile for contours, return the contours
% in the ISIdata struct.
% If the saveFig flag is set, a .fig file is written out with just the
% smoothed signal image

if ~isfield(ISIdata,'vesselmask')
    warning('isolatebarrel:nomask','No vessel mask in dataset.');
end
if nargin < 4 || isempty(saveFig) %selectROI flag doesn't exist
    saveFig = 0;
end

%% load and smooth signal frame

%select analysis ROI if specified
if isfield(ISIdata,'analysisROI')
    signalFrame = signalFrame.*(ISIdata.analysisROI.ROI);
end

%use vesselmask, set to ones if not using
filtframe=signalFrame.*ISIdata.vesselmask;

%crop data to avoid edge artifacts
croppix=0;
if croppix>0
    dipFramesAvg=filtframe(croppix:(ISIdata.frameSizeYX(2)-croppix),croppix:(ISIdata.frameSizeYX(1)-croppix));
else dipFramesAvg=filtframe; 
end

hfig = figure('name',['BARREL FINDER- ' prmtsCurrent.name]);

% F = fspecial('gaussian', 6, 2); %6x6 gaussian kernel, sigma=2
F = fspecial('gaussian', 12, 4); %12x12 gaussian kernel, sigma=4

filtframeSmooth1 = filter2(F, dipFramesAvg,'valid');
dipIm=filtframeSmooth1;
% filtframeSmooth2 = filter2(F, filtframeSmooth1,'valid');
% dipIm=filtframeSmooth2;

ISIdata.signalFrameFilt=dipIm;

himg=imagesc(dipIm);     
axis image; hold on; colormap hot; colorbar
title(prmtsCurrent.Whisker,'interpreter','none');
caxis([-6e-4 6e-4]);

if saveFig
    savefilename=fullfile(prmtsCurrent.path2dir,prmtsCurrent.name);
    [p n e]=fileparts(savefilename);
    
    %saveas(gcf,fullfile(p,[n '.fig']),'fig');
    saveas(gcf,fullfile(p,[n '.png']),'png');
end

%% get contours

contour_level = prmtsCurrent.ContourLineVals;
contour_level = sort(contour_level,'ascend');
nContourLevels = numel(contour_level);

if nContourLevels>1
    caxis(contour_level([1 end])-1)
else
    caxis([-4e-4 4e-4]);
end

% colorrange=[.00010 .00015 .0002 .00025 .0003 .00035];
contour_color={'r' 'g' 'b' 'c' 'm' 'y' 'k' 'r+-' 'g+-' 'b+-'};
nContColor=length(contour_color);

%ensure contour levels are in ascending order (lower value -> deoxigenated area)

contourLinesLeftPerLevel = zeros(nContourLevels,1);
gotLowest = 0;
for iLevel=1:nContourLevels
    [xcontour{1,iLevel},h]=contour(dipIm,[-1 -1]+contour_level(iLevel),contour_color{mod(iLevel,nContColor)+1}, 'LineWidth',2);

    %%%%%delete small objects
    hChildren = get(h,'children');

    for iChild = 1 : numel(hChildren)
        np = numel(get(hChildren(iChild),'Faces'));
        if np < prmtsCurrent.minFaces || np > prmtsCurrent.maxFaces
            delete(hChildren(iChild));
        end
    end

    %Check how many contour points left for current level
    hChildren = get(h,'children');
    contourLinesLeftPerLevel(iLevel) = numel(hChildren);

    %keep lowest for later use
    if contourLinesLeftPerLevel(iLevel) > 0
        if ~gotLowest 
            XdataLowest = get(hChildren,'Xdata');
            YdataLowest = get(hChildren,'Ydata');
            if isnumeric(XdataLowest) %if only one contour, still display to give chance to skip
                XdataLowest={XdataLowest};
                YdataLowest={YdataLowest};
            end
            if iscell(XdataLowest) %only one contour
                %                 need to ask user to choose
                %                 keyboard
                %mjp 2011.07.29 added ability to select more than one contour
                idx = ISI_selectContour(XdataLowest,YdataLowest,prmtsCurrent,dipIm,contour_level(iLevel)-1);
                drawnow;
                if idx > 0
                    Xcontour = XdataLowest{idx(1)};
                    Ycontour =YdataLowest{idx(1)};
                    if numel(idx)>1
                        for ix=2:length(idx)
                            Xcontour=[Xcontour; NaN; XdataLowest{idx(ix)};];
                            Ycontour=[Ycontour; NaN; YdataLowest{idx(ix)};];
                        end
                        Xcontour=[Xcontour; NaN];
                        Ycontour=[Ycontour; NaN];
                    end
                    gotLowest = 1;
                    lowestLevel = iLevel;
                end
            else
                gotLowest = 1;
                lowestLevel = iLevel;
                Xcontour=XdataLowest;
                Ycontour=YdataLowest;
            end

        end
        %clabel(xcontour{1,iLevel},h)
    end %extracting lowest xy data

%     axis square
%     axis on;grid on;
%     %
%     axis ij

end
close(hfig);

hfig=figure('name',['BARREL FINDER- ' prmtsCurrent.name]);
himg=imagesc(dipIm);     
if nContourLevels>1
    caxis(contour_level([1 end])-1)
else
    caxis([-4e-4 4e-4]);
end
axis image; hold on; colormap hot; colorbar
title(prmtsCurrent.Whisker,'interpreter','none');

%% now display lowest contour
if isfield(prmtsCurrent,'refImage')
    if ~isempty(prmtsCurrent.refImage)
        figure('name',['BARREL LOCATOR REF- ' prmtsCurrent.name]  )
        imRef = imread(fullfile(prmtsCurrent.path2dir,prmtsCurrent.refImage));

        %Compute size difference between filtered "dip" image and reference
        [deltaRC] = size(imRef) - size(dipIm);
        deltaR = fix(deltaRC(1)/2);
        deltaC = fix(deltaRC(2)/2);
        
        imshow(imRef)
        axis image
        hold on
        x = Xcontour + deltaR;
        y = Ycontour + deltaC;
        ISIdata.contourM=[x y];
        
        plot(x,y,'r-','LineWidth',2)
        title([prmtsCurrent.name ' response level = ' num2str(contour_level(lowestLevel)-1)],'interpreter','none');
        mx=mean(x(~isnan(x)));
        my=mean(y(~isnan(y)));
        %mjp 2011.09.22 should be centering txt from mean; for now it's offset
        text(mx , my , prmtsCurrent.Whisker,'color','r','fontSize',10);
        figname=get(gcf,'name');
        saveas(gcf,fullfile(prmtsCurrent.path2dir, [figname(1:end-4) '.pdf']),'pdf');
    end %processing image
else
    fprintf('\nNo reference image for %s', prmtsCurrent.name);
end

return