% Perform affine transformation of blue IOS image into the coordinate
% system of the GalvoScanner blue/vessel image
function BlueImageRegistration

debug = 0;

% Load IOS blue image
[sFile, sPath] = uigetfile('*.png', 'Select IOS vessel image');
sIOSBlueFile = [sPath sFile];
unregistered = imread(sIOSBlueFile);
unregistered = unregistered(:,:,1);

% Load the GalvoScanner blue/vessel image WITHOUT stereotactic coordinates
[sFile, sPath] = uigetfile('*.png', 'Select GalvoScanner vessel image wo/coordinates');
sGalvoScanBlueFile = [sPath sFile];
reference = imread(sGalvoScanBlueFile);
%reference = reference(:,:,1);

% Load the GalvoScanner blue/vessel image WITH stereotactic coordinates
[sFile, sPath] = uigetfile('*.png', 'Select GalvoScanner vessel image w/coordinates');
sGalvoScanBlueFileWCoords = [sPath sFile];
coordsimg = imread(sGalvoScanBlueFileWCoords);

% Load barrel map
[sFile, sPath] = uigetfile('*.fig', 'Select the barrel map (in same coordinates as IOS image)');
sBarrelMapFig = [sPath sFile];
hBarrelmap = open(sBarrelMapFig);
barrelmap = get(findobj(hBarrelmap, 'type', 'image'), 'cdata');
%barrelmap = barrelmap .* -1;

% Threshold and normalize barrel map (0-255);
barrelmap(barrelmap(:) < -1) = -1;
barrelmap(barrelmap(:) > -.5) = -.5;
barrelmap = barrelmap - min(barrelmap(:));
barrelmap = barrelmap ./ max(barrelmap(:));

% Check if controlpoint data already exists on disk
if exist('ImageReControlPoints.mat', 'file')
    sAns = questdlg('Do you want to load existing controlpoint data from disk?', ...
        'Load controlpoints', 'Yes', 'No', 'Cancel', 'Yes');
    switch sAns
        case 'Yes',
            load('ImageReControlPoints.mat', '-MAT');
        case 'Cancel',
            return;
    end
end

% Get control points interactively
if ~exist('input_points', 'var')
    [input_points, base_points] = cpselect(unregistered, reference, 'Wait', true);
    % Find transform
    TFORM = cp2tform(input_points, base_points, 'affine');
    % Save controlpoints and transform
    save ImageReControlPoints.mat input_points base_points TFORM
end

% Perform affine transformation of unregistered image (IOS blue/vessel image)
[nrows ncols] = size(reference);
registered = imtransform(unregistered, TFORM, 'XData',[1 ncols], 'YData',[1 nrows]);
registered = imresize(registered, [nrows ncols]);

% Transform barrel map
barrelmap_reg = imtransform(barrelmap, TFORM, 'XData',[1 ncols], 'YData',[1 nrows], 'FillValues', 1);
barrelmap_reg = imresize(barrelmap_reg, [nrows ncols]);

% Show all images
if debug
    figure
    subplot(1,3,1)
    imagesc(unregistered)
    axis image
    
    subplot(1,3,2)
    imagesc(reference)
    axis image
    
    subplot(1,3,3)
    imagesc(registered)
    axis image
end

% Click and mark stereotactic coordinates (should be in square grid, so two
% points is enough).
hFig = figure;
imagesc(coordsimg)
% Click on [-1 -3]
title('Click on (-1,-3)')
[nX1, nY1] = ginput(1);
% Click on [-2 -3]
title('Click on (-2,-3)')
[nX2, nY2] = ginput(1);
close(hFig)


% Show barrel map superimposed on galvoscanner vessel image
figure
%h1 = imshow(reference);
h1 = imshow(coordsimg);

hold on
h2 = imshow(barrelmap_reg);

alphamap = zeros(size(barrelmap_reg));
alphamap(barrelmap_reg < .8) = .3;
set(h2, 'alphadata', alphamap, 'alphaDataMapping','none');
colormap gray



return