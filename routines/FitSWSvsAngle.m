
%==========================================================================

function out = FitSWSvsAngle(cvalsToFit,anglesToFit,fixedphi,fixedthetatilt,setupSWS)

% Fit input SWS values to SWS ellipse expected in a TI material.
% cvalsToFit contains SWS in each direction
% anglesToFit contains the acquired rotation angles (corresponding to cvals)
% fixedphi is either a string labeled 'vary' or a double equal to the angle at which the phi should be fixed

% fixedthetatilt is typically set to 0 (no tilt correction), but can also
% be another number or a structure, to correct for tilt in various ways.
% (see seperate script, CvalsFromSqrtEllipse_IncTilt)

% setupSWS is a structure containing instructions for fitting
% out is a structure containing results of fit

if size(cvalsToFit,1)>1
    if size(cvalsToFit,2)==1
        cvalsToFit=cvalsToFit'; % flip if oriented incorrectly-- needs to be 1x##
    end
end

cParvals  = (setupSWS.minspeed*2):0.1:setupSWS.maxspeed;      % values to try "parallel" to fibers
cPerpvals = (setupSWS.minspeed):0.1:(setupSWS.maxspeed*.5);   % values to try "perpendicular" to fibers

if or(isstring(fixedphi),ischar(fixedphi))
    if strcmp(fixedphi,'vary')
        phivals = -20:1:20; % values of ellipse rotation to try
    elseif strcmp(fixedphi,'vary90to90')
        phivals=-90:1:89;
    elseif regexp(fixedphi,'varyaround.*')
        tmp=regexp(fixedphi,'varyaround(.*)','tokens');
        centerphi=str2double(tmp{1}{1});
        phivals=centerphi+(-10:1:10);
    end
else
    phivals=fixedphi; %fix at the previously determined rotation value
end
if or(isstring(fixedthetatilt),ischar(fixedthetatilt))
    error('moving tilt angle not incoporated')
else
    thetatiltvals=fixedthetatilt; %pass through for tilt values
end

cParRange  = max(cParvals) -min(cParvals);
cPerpRange = max(cPerpvals)-min(cPerpvals);
phiRange   = max(phivals)  -min(phivals);

% find an initial "best guess" from course spacing of params to try
[cParfit,cPerpfit,phifit,~] = FitSWSvsAngleVaryParams(cvalsToFit,anglesToFit,cParvals,cPerpvals,phivals,thetatiltvals);

% allow for additional movement in phi, though not in other two fit
% parameters
if or(isstring(fixedphi),ischar(fixedphi)) %onlyshift if fixedphi="vary", not if number
    for ishiftphi=1:8 % 8 lets get to +/- 100 at one extreme, if somehow at rotation of 90
        if phifit==phivals(1)
            phivals=phivals-10; % shift by -10 degrees
        elseif phifit==phivals(end)
            phivals=phivals+10; % shift by 10 degrees
        else
            break % phifit not at edge of tested values, so move on
        end
        [cParfit,cPerpfit,phifit,~] = FitSWSvsAngleVaryParams(cvalsToFit,anglesToFit,cParvals,cPerpvals,phivals,thetatiltvals);
    end
end


niterations=5;
for ii=1:niterations            % loop to iteratively refine search range, narrow in on best params

    cParRange  = cParRange /2;
    cPerpRange = cPerpRange/2;
    phiRange   = phiRange  /2;

    [cParvals,cPerpvals,phivals] = RefineParameterGuesses(cParfit,cPerpfit,phifit,cParRange,cPerpRange,phiRange);

    [cParfit,cPerpfit,phifit,cfit] = FitSWSvsAngleVaryParams(cvalsToFit,anglesToFit,cParvals,cPerpvals,phivals,thetatiltvals);
end

%output cvalues of entire range fit
anglestocompute=-90:270;
cfitfull = CvalsFromSqrtEllipse_IncTilt(anglestocompute,cParfit,cPerpfit,phifit,thetatiltvals);

out.cvalsToFit   = cvalsToFit;
out.anglesToFit  = anglesToFit;
out.cParfit      = cParfit;
out.cPerpfit     = cPerpfit;
out.phifit       = phifit;
out.thetatilt    = fixedthetatilt;
out.cfit         = cfit;
out.anglesfull   = anglestocompute;
out.cfitfull     = cfitfull;

end

%==========================================================================

function [cParfit,cPerpfit,phifit,cfit] = FitSWSvsAngleVaryParams(cvalsToFit,anglesToFit,cParvals,cPerpvals,phivals,thetatiltvals)

ssd = NaN(length(cParvals),length(cPerpvals),length(phivals));

cvalsToFitmatrix=permute(repmat(cvalsToFit,length(cParvals),1,length(cPerpvals)),[1 3 2]);     % replicate, angle is third dimension

if length(thetatiltvals)>1
    error('tilt value should be set, many not coded')
end

for iphi = 1:length(phivals)
    %try each combo of cpar, cperp, iphi

    cvals = CvalsFromSqrtEllipse_IncTilt_Matrix(anglesToFit,cParvals,cPerpvals,phivals(iphi),thetatiltvals); %cvals for that combo of params
    ssd(:,:,iphi) = sum(((cvalsToFitmatrix-cvals).^2),3,'omitnan'); % sum squared difference between guess and data-- sum across cvals, 3rd dimension
end

cParvalsmatrix=repmat(cParvals',1,length(cPerpvals),length(phivals));
cPerpvalsmatrix=repmat(cPerpvals,length(cParvals),1,length(phivals));
ssd(cPerpvalsmatrix>cParvalsmatrix)=NaN; % exclude Cperp>cPar, likely a bad phi fit

[~,minidx] = min(ssd(:)); %look for minimum
[icParmin,icPerpmin,iphimin]=ind2sub(size(ssd),minidx); %find which index is minimal
cParfit = cParvals(icParmin); % pick the params corresponding to the minimal fit
cPerpfit = cPerpvals(icPerpmin);
phifit = phivals(iphimin);

cfit = CvalsFromSqrtEllipse_IncTilt(anglesToFit,cParfit,cPerpfit,phifit,thetatiltvals); %report out the best case cvals
end

%==========================================================================

function [newpar1vals,newpar2vals,newpar3vals] = RefineParameterGuesses(par1,par2,par3,par1range,par2range,par3range)

% set new parameter search space as half of previous range, centered about "best fit" guess from previous iteration

par1min = par1 - par1range/2;
par1max = par1 + par1range/2;

par2min = par2 - par2range/2;
par2max = par2 + par2range/2;

par3min = par3 - par3range/2;
par3max = par3 + par3range/2;

npts=9;
newpar1vals = par1min:(par1range/npts):par1max;
newpar2vals = par2min:(par2range/npts):par2max;
newpar3vals = par3min:(par3range/npts):par3max;

if and(isempty(newpar3vals),par3range==0)
    newpar3vals=par3;
end
end

%==========================================================================
