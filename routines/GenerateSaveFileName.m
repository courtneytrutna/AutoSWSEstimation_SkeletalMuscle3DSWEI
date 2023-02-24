function [savefolder,SWSsettingsname]=GenerateSaveFileName(dataDir,setupdataprocessing)

% Generates filename for saving ellipse fits based on data directory (dataDir) and
% settigns held in setupdataprocessing.

spacetimeplanefolder=['ForSpaceTimePlanes' setupdataprocessing.infostring.PlaneInfo];
SWStypefolder=['SWSMethod' setupdataprocessing.SWSEstimationMethod];
SWSsettingsname=['SWSSettings' setupdataprocessing.SWSEstimationParams setupdataprocessing.setrandsampoutlierthres]; % for readability, include underscore as first character if interpstring isn't just ''


savefolder=[dataDir '/SWS3DEstimation/' spacetimeplanefolder '/' SWStypefolder '/'];
if ~exist(savefolder,'dir');mkdir(savefolder);end

end