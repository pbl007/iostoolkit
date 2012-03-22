%% ISI_analysis_multi script - add files names to be analyzed 
% modified from ISI_analysis_trialBased_PabloXXXXXXXX.m
% mjp 2011.07.04


clear
clc
close all


%% define analysis behavior
%NB: could also move these into specific file structs
prmts.saveMovie = 0;
prmts.saveToMat = 0;
prmts.analyzeEcoG = 0; %not used
prmts.useNoLight = 0;
prmts.useVesselMask = 0;
prmts.saveFig = 0; %write smoothed signal frame to .fig


%% define parameters and keep them in prmts structure

% basepath='D:\dklab\Mike IOS data\';
basepath='E:\IOS data\110111.Per.BS37';

subpath='';

path2dir=fullfile(basepath,subpath);
default = struct('path2dir',path2dir,'precision','int16','name',{[]},...
    'Whisker',{[]},'Trials2Use',[],'ECoGFile',[],'preStimDurSec',4,'stimDurSec',4,'Thresh',1e-4,...
    'ContourLineVals',0.9995:1e-4:1.0005,'minFaces',50,'maxFaces',1000,'selectROI',0);
%NB: ContourLineVals, minFaces, and maxFaces will all need to be changed to
%match the data set
default.ContourLineVals=0.9995:.5e-4:1.0005;

%analysis tips: 
%if not offering contour selection with the contour you
%  want, change the minFaces and maxFaces parameters to include the desired
%  contour(s). higher minFaces will show fewer small contours
%if want to change max and min on colorbar, change the contourLineVals min
%  and max

default.stimInterval=[0.5 5]; %default to 0.5-2 seconds of integration
default.maskthresh=0.45;


%% nolight should be calculated with the illumination light turned off and
% no whisker deflection
prmts.Rnolight=default;
prmts.Rnolight.name='';

%% vesselfile
vesselfile='blue.png';

%% file specific information
%Note: leave Trials2Use empty to use all

fqi = 0;%file counter

%insert file specific options
fqi = fqi + 1;
prmts.filesQueue(fqi) = default;

prmts.filesQueue(fqi).name = 'C1_153421.dat';   % !!!!
prmts.filesQueue(fqi).Whisker = 'C1';           % !!!!

prmts.filesQueue(fqi).refImage = vesselfile;
prmts.filesQueue(fqi).selectROI =0;
prmts.filesQueue(fqi).Trials2Use =[];

%prmts.filesQueue(fqi).contourQuantile=0.05;


% fqi = fqi + 1;
% prmts.filesQueue(fqi) = default;
% prmts.filesQueue(fqi).name = '20090618-B4';
% prmts.filesQueue(fqi).Whisker = 'B4';
% prmts.filesQueue(fqi).ECoGFile = '';
% prmts.filesQueue(fqi).refImage = 'Blue_02.png';
% prmts.filesQueue(fqi).preStimDurSec = 1;
% prmts.filesQueue(fqi).stimDurSec = 1;
% prmts.filesQueue(fqi).Trials2Use =[];

%% run 

prmts = ISI_analysis(prmts);



