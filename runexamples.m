% runexamples.m

% Code that runs scripts over example files in git repo:
% https://github.com/courtneytrutna/AutoSWSEstimation_SkeletalMuscle3DSWEI. 

% assumes particular directory structure and format of data loaded for
% analysis, see README on gitrepo for more details. 

% Tested in MATLAB 2022a. Various additional MATLAB Toolboxes may be required.

%% Setup

% List data directories
dataDir=dir('SampleData/InVivo/'); % this directory should be where your data is stored
%dataDir=dir('SampleData/Phantom/'); % this directory should be where your data is stored
dataDir = dataDir(~startsWith({dataDir.name},'.')); % remove hidden files

% Information about your acquistion
setupdataprocessing.infostring.PlaneInfo=''; %use to distingush the same data preprocessed in different ways
setupdataprocessing.anglesDeg=-90:5:265; % what rotation angles are in your dataset?

% Set fitting method and lateral range
SWSMethodFlag=1; % values of 1 to 15 depending on desired setting, see setSWSfitparams function below
lateralrange=[2 18]; % set lateral range in mm. 2-18mm as default
%lateralrange=[4 18 2 16]; % use four element vector for angle informed methods, where second set is values across fibers
setupdataprocessing=setSWSfitparams(setupdataprocessing,SWSMethodFlag,lateralrange); % additional parameters for fitting. see below

%% Run Fits
overwriteexistingflag=1; % set to zero to skip datasets that have already been processed
savelightweightflag=0; % set to zero to only save key outputs as .txt file, otherwise .mat and figure created
for i=1:length(dataDir)
    disp(['Running Over Data: ' dataDir(i).name])
    RunAuto3DSWSFromPlanes([dataDir(i).folder '/' dataDir(i).name],setupdataprocessing,overwriteexistingflag,savelightweightflag)
end

%==========================================================================

%% Function for Setting Fit Parameters

function setupdataprocessing=setSWSfitparams(setupdataprocessing,SWSMethodFlag,lateralrange)

latmin=lateralrange(1);
latmax=lateralrange(2);
if length(lateralrange)==4
    latmin_across=lateralrange(3);
    latmax_across=lateralrange(4);
end

% Radon/TTP/XCorr parameter settings determined to be optimium
if any(SWSMethodFlag==1:5) % Radon
    originthres=2;
    multiwavedetectionflag='_multiwave';
    %multiwavedetectionflag=''; % to skip multiwave detection. Note For Radon Sum, will still find, but will not analyze
elseif any(SWSMethodFlag==6:10) % TTP
    originthres=1;
    ttpoutlier=1;
    multiwavedetectionflag='_multiwave';
    %multiwavedetectionflag=''; % to skip multiwave detection
elseif any(SWSMethodFlag==11:15) % XCorr
    xcorrsigthres=0;
    xcorrlatstepmm=4;
    xcorrcoeffthres=0.5;
    xcorroutlier=12;
end

% Dyanmic fits settings (determined empirically)
dynminlat=3;
ttpr2thres=0.9;
xcorrdynerrorthres=0.1;

switch SWSMethodFlag
    case 1
        setupdataprocessing.SWSEstimationMethod='_RadonSumFixedLatRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_lat' num2str(latmin) 'to' num2str(latmax) 'originthres' num2str(originthres) multiwavedetectionflag];
    case 2
        setupdataprocessing.SWSEstimationMethod='_RadonSumAngleInformedLatRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_Alonglat' num2str(latmin) 'to' num2str(latmax) 'Acrosslat' num2str(latmin_across) 'to' num2str(latmax_across) 'originthres' num2str(originthres) multiwavedetectionflag];
        setupdataprocessing.fiberestimatefile='/FiberAngleEstimate/fiber_orientation.mat';
    case 3
        setupdataprocessing.SWSEstimationMethod='_RadonSumDynamicLatRangeRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_originthres' num2str(originthres) 'minlatrange' num2str(dynminlat)];
    case 4
        setupdataprocessing.SWSEstimationMethod='_RadonSumFixedLat3DFit';
        setupdataprocessing.SWSEstimationParams=['_lat' num2str(latmin) 'to' num2str(latmax) 'originthres' num2str(originthres)];
    case 5
        setupdataprocessing.SWSEstimationMethod='_RadonSumAngleInformedLat3DFit';
        setupdataprocessing.SWSEstimationParams=['_Alonglat' num2str(latmin) 'to' num2str(latmax) 'Acrosslat' num2str(latmin_across) 'to' num2str(latmax_across) 'originthres' num2str(originthres)];
        setupdataprocessing.fiberestimatefile='/FiberAngleEstimate/fiber_orientation.mat';
    case 6
        setupdataprocessing.SWSEstimationMethod='_TTPFixedLatRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_lat' num2str(latmin) 'to' num2str(latmax) 'originthres' num2str(originthres) 'outlierTOFthres' num2str(ttpoutlier) multiwavedetectionflag];
    case 7
        setupdataprocessing.SWSEstimationMethod='_TTPAngleInformedLatRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_Alonglat' num2str(latmin) 'to' num2str(latmax) 'Acrosslat' num2str(latmin_across) 'to' num2str(latmax_across) 'originthres' num2str(originthres) 'outlierTOFthres' num2str(ttpoutlier) multiwavedetectionflag];
    case 8
        setupdataprocessing.SWSEstimationMethod='_TTPDynamicLatRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_originthres' num2str(originthres) 'outlierTOFthres' num2str(ttpoutlier) 'minlatrange' num2str(dynminlat) 'R2thres' num2str(ttpr2thres)];
    case 9
        setupdataprocessing.SWSEstimationMethod='_TTPFixedLat3DFit';
        setupdataprocessing.SWSEstimationParams=['_lat' num2str(latmin) 'to' num2str(latmax) 'originthres' num2str(originthres) 'outlierTOFthres' num2str(ttpoutlier)];
    case 10
        setupdataprocessing.SWSEstimationMethod='_TTPAngleInformedLat3DFit';
        setupdataprocessing.SWSEstimationParams=['_Alonglat' num2str(latmin) 'to' num2str(latmax) 'Acrosslat' num2str(latmin_across) 'to' num2str(latmax_across) 'originthres' num2str(originthres) 'outlierTOFthres' num2str(ttpoutlier)];
    case 11
        setupdataprocessing.SWSEstimationMethod='_XCorrFixedLatRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_lat' num2str(latmin) 'to' num2str(latmax) 'outlierthres' num2str(xcorroutlier) 'Xcorrcoeffthres' num2str(xcorrcoeffthres) 'latstepmm' num2str(xcorrlatstepmm) 'sigthreshold' num2str(xcorrsigthres)];
    case 12
        setupdataprocessing.SWSEstimationMethod='_XCorrAngleInformedLatRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_Alonglat' num2str(latmin) 'to' num2str(latmax) 'Acrosslat' num2str(latmin_across) 'to' num2str(latmax_across) 'outlierthres' num2str(xcorroutlier) 'Xcorrcoeffthres' num2str(xcorrcoeffthres)  'latstepmm' num2str(xcorrlatstepmm) 'sigthreshold' num2str(xcorrsigthres)];
        setupdataprocessing.fiberestimatefile='/FiberAngleEstimate/fiber_orientation.mat';
    case 13
        setupdataprocessing.SWSEstimationMethod='_XCorrDynamicLatRANSACEllipse';
        setupdataprocessing.SWSEstimationParams=['_outlierthres' num2str(xcorroutlier) 'Xcorrcoeffthres' num2str(xcorrcoeffthres)  'latstepmm' num2str(xcorrlatstepmm) 'sigthreshold' num2str(xcorrsigthres) 'minlatrange' num2str(dynminlat) 'Xcorrmeanerrorthres' num2str(xcorrdynerrorthres)];
    case 14
        setupdataprocessing.SWSEstimationMethod='_XCorrFixedLat3DFit';
        setupdataprocessing.SWSEstimationParams=['_lat' num2str(latmin) 'to' num2str(latmax) 'outlierthres' num2str(xcorroutlier) 'Xcorrcoeffthres' num2str(xcorrcoeffthres)  'latstepmm' num2str(xcorrlatstepmm) 'sigthreshold' num2str(xcorrsigthres)];
    case 15
        setupdataprocessing.SWSEstimationMethod='_XCorrAngleInformedLat3DFit';
        setupdataprocessing.SWSEstimationParams=['_Alonglat' num2str(latmin) 'to' num2str(latmax) 'Acrosslat' num2str(latmin_across) 'to' num2str(latmax_across) 'outlierthres' num2str(xcorroutlier) 'Xcorrcoeffthres' num2str(xcorrcoeffthres)  'latstepmm' num2str(xcorrlatstepmm) 'sigthreshold' num2str(xcorrsigthres)];
        setupdataprocessing.fiberestimatefile='/FiberAngleEstimate/fiber_orientation.mat';
end

if regexp(setupdataprocessing.SWSEstimationMethod,'RANSACEllipse')
    setupdataprocessing.setrandsampoutlierthres='';%'_randsampoutthres0.5'; % leave plane for default value of 1
end

end
