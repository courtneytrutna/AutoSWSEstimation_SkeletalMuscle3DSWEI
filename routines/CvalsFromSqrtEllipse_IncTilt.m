function cvals = CvalsFromSqrtEllipse_IncTilt(anglesDeg,cPar,cPerp,phi,tiltAngleDeginput)
% calculates c values of ellipse (as a function of anglesDeg) for the
% given parameters (cPar,CPerp,phi,tiltAngleDeginput)

% note we typically do not assume the fibers are substanitally tilted out
% of the plane, and use tiltAngleDeginput as 0. Using tilt correction can
% introduce another source of error if tilt estimation is not accurate.

% if the fibers are tilted, there are two potential ways to correct:
% modifying the angles of the ellipse, or using a cosine correction only in
% the major axis direction. Which is more appropraite depends on the height
% of the push and if elliptical-like propogation is visualized in the
% axial-lateral plane (which would call for modifying the ellipse angles)
% or if the pattern is more like a cylindrically spreading wave (calling
% for a cosine correction).

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
anglesReAxis = acosd(cosd(tiltAngleDeg)*cosd(anglesDeg-phi));


muVsTheta = c44*c66 ./ ( c44*(sind(anglesReAxis)).^2 + c66*(cosd(anglesReAxis)).^2 );
cvals = sqrt(muVsTheta/rho);
end
