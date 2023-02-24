
%==========================================================================

function [out,removedoutlierflag] = FitSWSvsAngle_RANSAC(cvalsToFit,anglesToFit,fixedphi,fixedthetatilt,fitparams,setupSWS)

% Performs RANSAC fitting of ellipse fit, calling FitSWSvsAngle many times
% with random subselection of datapoints and calculating cost function

% cvalsToFit contains SWS in each direction
% anglesToFit contains the acquired rotation angles (corresponding to cvals)
% fixedphi is either a string labeled 'vary' or a double equal to the angle at which the phi should be fixed

% fixedthetatilt is typically set to 0 (no tilt correction), but can also
% be another number or a structure, to correct for tilt in various ways.
% (see seperate script, CvalsFromSqrtEllipse_IncTilt)

% fitparams is a strucutre containing instructions for fitting
% setupSWS is a structure containing instructions for fitting
% out is a structure containing results of fit
% removedoutlierflag contains 0's and 1's denoting angles that were detected as outliers (1's)

if all(isnan(cvalsToFit))
    out.cvalsToFit   = cvalsToFit;
    out.anglesToFit  = anglesToFit;
    out.cParfit      = NaN; % give nan because failed fit
    out.cPerpfit     = NaN;
    out.phifit       = NaN;
    out.thetatilt    = fixedthetatilt;
    out.cfit         = NaN;
    out.anglesfull   = NaN;
    out.cfitfull     = NaN;
    out.ellipsecostfunctionval=NaN;
    out.ellipsecost_inlier=NaN;
    removedoutlierflag=ones(size(cvalsToFit));

    return
end

niter=100;

currentbestellipsequal=length(cvalsToFit)*10; % this is a terrible error, and will almost certainly be passed.

rng('default'); % reset for reproducibility
for iiterations=1:niter
    itest=randi(length(cvalsToFit),5,1);

    outFitSWSvsAngleSH=FitSWSvsAngle(cvalsToFit(itest),anglesToFit(itest),fixedphi,fixedthetatilt,setupSWS);
    outFitSWSvsAngleSH.cvalsToFit=cvalsToFit;outFitSWSvsAngleSH.anglesToFit=anglesToFit;

    % remove outliers:
    [~,~,weightedellipsequal_outlieratthres,weightedellipsequal_inliermeansqerr]=RemoveOutliers_CalcErrorEllipse(outFitSWSvsAngleSH,fitparams.randsampoutlierthreshold,fixedthetatilt);

    if weightedellipsequal_outlieratthres<currentbestellipsequal % use RANSAC for fitting
        currentbestellipsequal=weightedellipsequal_outlieratthres;
        currentbestellipsequal_inlier=weightedellipsequal_inliermeansqerr; %carry inlier only for later quality metric imposition
        currentbestout=outFitSWSvsAngleSH;
    end

end

[speed_bestfit_cleaned,removedoutlierflag,~,~]=RemoveOutliers_CalcErrorEllipse(currentbestout,fitparams.randsampoutlierthreshold,fixedthetatilt);

out=FitSWSvsAngle(speed_bestfit_cleaned,anglesToFit,fixedphi,fixedthetatilt,setupSWS);

out.ellipsecostfunctionval=currentbestellipsequal;
out.ellipsecost_inlier=currentbestellipsequal_inlier;
end