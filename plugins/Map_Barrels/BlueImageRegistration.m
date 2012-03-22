% Perform affine transformation of blue IOS image into the coordinate
% system of the GalvoScanner blue/vessel image
transformtype  = 'affine';
unregistered = imread('IOS_Blue.png');
unregistered = unregistered(:,:,1);
reference = imread('GalvoScan_Blue_LoContrast.png');
%reference = reference(:,:,1);

% Load barrel map
hBarrelmap = open('BS41_BarrelMap_ROINormalized_Filtered.fig');
barrelmap = get(findobj(hBarrelmap, 'type', 'image'), 'cdata');
%barrelmap = barrelmap .* -1;

% Threshold and normalize barrel map (0-255);
barrelmap(barrelmap(:) < -1) = -1;
barrelmap(barrelmap(:) > -.5) = -.5;
barrelmap = barrelmap - min(barrelmap(:));
barrelmap = barrelmap ./ max(barrelmap(:));

[input_points, base_points] = cpselect(unregistered, registered, 'Wait', true);

% Save corresponding points
save ImageReControlPoints.mat input_points base_points TFORM

% Find transform 
TFORM = cp2tform(input_points, base_points, transformtype);
 
% Perform affine transformation of unregistered image (IOS blue/vessel image)
[nrows ncols] = size(reference);
registered = imtransform(unregistered, TFORM, 'XData',[1 ncols], 'YData',[1 nrows]);
registered = imresize(registered, [nrows ncols]);

% Transform barrel map
barrelmap_reg = imtransform(barrelmap, TFORM, 'XData',[1 ncols], 'YData',[1 nrows], 'FillValues', 1);
barrelmap_reg = imresize(barrelmap_reg, [nrows ncols]);

% Show all images
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

% Show barrel map on top of galvoscanner vessel image
figure
h1 = imshow(reference);
hold on
h2 = imshow(barrelmap_reg);
alphamap = zeros(size(barrelmap_reg));
alphamap(barrelmap_reg < .8) = .3;
set(h2, 'alphadata', alphamap, 'alphaDataMapping','none');
colormap jet

