function [cvals_cleaned,removeflag,ellipsequal_outlieratthres,ellipsequal_inliermeansqerr]=RemoveOutliers_CalcErrorEllipse(out_ellipsefit,diffmax,fixedthetatilt)
% determine RANSAC cost function by first detecting large errors, then
% summing number of outliers + square error of inliers
% out_ellipsefit is a strucutre contaiing information about the ellipse fit
% diffmax is the outlier threshold value (in terms of error)
% fixedthetatilt contains information about fitting (need to pass through because structure not maintained in out_ellipse)

% cvals_cleaned is the output SWS as a function of angle with outliers
% removeflag is a vector the same size as cvals_cleaned, 1 for outliers, 0 for inliers
% two quality metrics are output:
% ellipsequal_outlieratthres, the RANSAC cost function and
% ellipsequal_inliermeansqerr, used for overal data quality exclusion

cvals = CvalsFromSqrtEllipse_IncTilt(out_ellipsefit.anglesToFit,out_ellipsefit.cParfit,out_ellipsefit.cPerpfit,out_ellipsefit.phifit,fixedthetatilt); %cvals for that combo of params
difference = ((out_ellipsefit.cvalsToFit-cvals).^2); % sum squared difference between guess and data-- sum across cvals, 3rd dimension

removeflag=difference>diffmax;
cvals_cleaned=out_ellipsefit.cvalsToFit;
cvals_cleaned(removeflag)=NaN;

% output a ransac cost function:
% C = sum across all cval({weighteddiff if weighteddiff<weightdiffmax; weightdiffmax is weighteddiff>weighteddiff})
% See Wang et al, UMB 2010

difference(removeflag)=diffmax;
ellipsequal_outlieratthres=sum(difference,'omitnan');

% recalc without the outliers and return both
difference(removeflag)=NaN;
ellipsequal_inliermeansqerr=mean(difference,'omitnan');

end
