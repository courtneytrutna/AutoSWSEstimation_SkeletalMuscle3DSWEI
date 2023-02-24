function [latstart,latend,rot_angle]=setAngleInformedRange(latstartAlong,latendAlong,latstartAcross,latendAcross,setupdataprocessing,dataDir)
% loads information about fiber rotation angle and uses to set AngleInformed lateral range

% latstartAlong, latendAlong,latstartAcross, and latendAcross set ranges in
% the two principle directions. All other angles are assigned based off the
% ellipses using these start/end values as major and minor axes.

% dataDir and setupdataprocessing together dictate the location and name of the fiber estimate file to use
% this file should have a variable named rot_angles, with a single number
% representing the angle of fiber rotation relative to the 0 degrees angle
% of setupdataprocessing.anglesDeg

% latstart and latend contain start and end values for each angle
% rot_angle contains the fiber rotation angle

fiberangleest=load([dataDir setupdataprocessing.fiberestimatefile]);

%ellipse equation is r=a*b/sqrt((b*cosd(theta))^2+(a*sind(theta))^2)
thetarange=setupdataprocessing.anglesDeg-fiberangleest.rot_angles;
rot_angle=fiberangleest.rot_angles;

latstart=latstartAlong*latstartAcross./sqrt((latstartAcross*cosd(thetarange)).^2+(latstartAlong*sind(thetarange)).^2);
latend=latendAlong*latendAcross./sqrt((latendAcross*cosd(thetarange)).^2+(latendAlong*sind(thetarange)).^2);
end