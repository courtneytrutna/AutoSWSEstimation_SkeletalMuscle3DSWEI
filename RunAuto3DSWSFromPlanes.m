function RunAuto3DSWSFromPlanes(dataDir,setupdataprocessing,overwriteexistingflag,savelightweightflag)
% Main function for SWS estimation, this script:
% --> Loads data (space-time planes),
% --> runs appropriate scripts for SWS estimation based on input 
% --> saves SWS information in either lightweight or full form 

% dataDir gives location of file to run
% setupdataprocessing is a structure with many expected fields, see runexamples.m
% overwriteexistingflag=1 will skip fits that have already been saved, set to 0 to rerun and overwrite
% savelightweightflag=1 will save a small .txt file with select output parameters. Set to 0 to save additional information and plot fit results

if nargin<3
    overwriteexistingflag=0;
end
if nargin<4
    savelightweightflag=0;
end

close all;
addpath routines 

[savefolder,SWSsettingsname]=GenerateSaveFileName(dataDir,setupdataprocessing);

if savelightweightflag
    savefilename=[savefolder 'LightOverall3DSWS_' SWSsettingsname '.txt'];
else
    savefilename=[savefolder 'Overall3DSWS_' SWSsettingsname '.mat'];
end

disp(['Running Method: ' setupdataprocessing.SWSEstimationMethod(2:end) SWSsettingsname])

if and(~overwriteexistingflag,exist(savefilename,'file'))
    % skip if already exists & don't want to overwrite
    disp('this file already exists, continuing on')
    return
end
%% find & load space time plane file
planeDir=[dataDir '/PlanesForAnalysis/'];
setup3DSWS.template = ['dataplanes' setupdataprocessing.infostring.PlaneInfo '.mat'];
filetoload= dir([planeDir setup3DSWS.template]);

tmp=load([filetoload.folder '/' filetoload.name]);
outplanes=tmp.out;
if isfield(tmp,'setupMakePlanes')
setupMakePlanes=tmp.setupMakePlanes;
else
    setupMakePlanes.climsforplane=[-1 1]; % set for plotting if information not saved with planes
    if ~savelightweightflag;disp('plotting with default colorscale');end
end

setupdataprocessing.nangles=length(setupdataprocessing.anglesDeg); % number of rotation angles in the dataset

%% Additional Global Settings
% other settings are implemented through setupdataprocessing or within the
% method-specific function

setup3DSWS.fignum=1*(savelightweightflag==0);
setup3DSWS.maxspeed=10;
setup3DSWS.minspeed=0.25;

%% Find SWS using set 3DSWEI algorithm

if setup3DSWS.fignum
% setup figures
figure(1);clf;tmp=gcf;tmp.Position=[1 1 1700 860];
end

switch setupdataprocessing.SWSEstimationMethod
    % Radon Sum Based Options
    case {'_RadonSumFixedLatRANSACEllipse','_RadonSumAngleInformedLatRANSACEllipse'}
        [outSWS,setup3DSWS]=FindSWS_RadonSum_SetLatRange_Each(outplanes,setupdataprocessing,setup3DSWS,setupMakePlanes,dataDir);
        [out,setup3DSWS]=FindSWS_EllipseFit(dataDir,setupdataprocessing,outSWS,setup3DSWS);
    case '_RadonSumDynamicLatRangeRANSACEllipse'
        [outSWS,setup3DSWS]=FindSWS_RadonSum_DynamicLatRange_Each(outplanes,setupdataprocessing,setup3DSWS,setupMakePlanes);
        [out,setup3DSWS]=FindSWS_EllipseFit(dataDir,setupdataprocessing,outSWS,setup3DSWS);
    case {'_RadonSumFixedLat3DFit','_RadonSumAngleInformedLat3DFit'}
        [out,setup3DSWS]=FindSWS_RadonSum_SetLatRange_3D(outplanes,setupdataprocessing,setup3DSWS,setupMakePlanes,dataDir);
       
        % TTP Based Options
    case {'_TTPFixedLatRANSACEllipse','_TTPAngleInformedLatRANSACEllipse','_TTPDynamicLatRANSACEllipse'}
        if regexp(setupdataprocessing.SWSEstimationParams,'.*multiwave.*')
            [outTOF,setup3DSWS]=Calc_TTP_Multiwave(outplanes,setup3DSWS);
            [outSWS,setup3DSWS]=FindSWS_TTP_Multiwave_Each(outTOF,outplanes,setupdataprocessing,setup3DSWS,setupMakePlanes,dataDir);
        else
            [outTOF,setup3DSWS]=Calc_TTP(outplanes,setup3DSWS);
            [outSWS,setup3DSWS]=FindSWS_TTP_Each(outTOF,outplanes,setupdataprocessing,setup3DSWS,setupMakePlanes,dataDir);
        end
        [out,setup3DSWS]=FindSWS_EllipseFit(dataDir,setupdataprocessing,outSWS,setup3DSWS);
    case {'_TTPFixedLat3DFit','_TTPAngleInformedLat3DFit'}
        [outTOF,setup3DSWS]=Calc_TTP(outplanes,setup3DSWS);
        [out,setup3DSWS]=FindSWS_TTP_SetLatRange_3D(outTOF,outplanes,setupdataprocessing,setup3DSWS,setupMakePlanes,dataDir);

        % Xcorr Based Options
    case {'_XCorrFixedLatRANSACEllipse','_XCorrAngleInformedLatRANSACEllipse','_XCorrDynamicLatRANSACEllipse'}
        [outXCorr,setup3DSWS]=Calc_XCorr(outplanes,setup3DSWS,setupdataprocessing);
        [outSWS,setup3DSWS]=FindSWS_XCorr_Each(outXCorr,outplanes,setupdataprocessing,setup3DSWS,setupMakePlanes,dataDir);
        [out,setup3DSWS]=FindSWS_EllipseFit(dataDir,setupdataprocessing,outSWS,setup3DSWS);
    case {'_XCorrFixedLat3DFit','_XCorrAngleInformedLat3DFit'}
        [outXCorr,setup3DSWS]=Calc_XCorr(outplanes,setup3DSWS,setupdataprocessing);
        [out,setup3DSWS]=FindSWS_XCorr_SetLatRange_3D(outXCorr,outplanes,setupdataprocessing,setup3DSWS,setupMakePlanes,dataDir);
    
end


%% Save File

if savelightweightflag
    if isfield(out,'costfunctionval') % only exists for EllipseFits
        writetable(table(out.cPar,out.cPerp,out.costfunctionval,out.costfunctionval_inlier,out.percentpts,'VariableNames',{'cPar','cPerp','EllipseCostVal','InlierMeanSqErr','PercentPtsIncluded'}),savefilename)
    else % for 3D-Fits, can't save
        writetable(table(out.cPar,out.cPerp,'VariableNames',{'cPar','cPerp'}),savefilename)
    end
    disp('lightweight save complete')
else
    save(savefilename,'out','setup3DSWS','-v7.3')
    disp('full save complete')
end

end
