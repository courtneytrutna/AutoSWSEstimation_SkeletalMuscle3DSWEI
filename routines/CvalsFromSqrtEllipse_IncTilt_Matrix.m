function cvals = CvalsFromSqrtEllipse_IncTilt_Matrix(anglesDeg,cPar,cPerp,phi,tiltAngleDeginput)
% calculates c values of ellipse (as a function of anglesDeg) for the
% given parameters (cPar,CPerp,phi)
% note coded to allow for cPar, cPerp to be matrixes, while other
% implementations assume single numbers for these values

if isstruct(tiltAngleDeginput)
    tiltcorrectionforCL=tiltAngleDeginput.tiltcorrectionforCL;
    tiltAngleDeg=tiltAngleDeginput.tiltforellipse; %fix at the previously determined tilt value
else
    tiltcorrectionforCL=0;
    tiltAngleDeg=tiltAngleDeginput;
end
cPar=cosd(tiltcorrectionforCL)*cPar;

rho = 1000;
c44 = rho * cPar.^2;
c66 = rho * cPerp.^2;
c44=repmat(c44',[1,length(cPerp)]);
c66=repmat(c66,[length(cPar),1]);

for i=1:length(anglesDeg)
    anglesReAxis = acosd(cosd(tiltAngleDeg)*cosd(anglesDeg(i)-phi));
    muVsTheta(:,:,i) = c44.*c66 ./ ( c44*(sind(anglesReAxis)).^2 + c66*(cosd(anglesReAxis)).^2 );
end
cvals = sqrt(muVsTheta/rho);
end