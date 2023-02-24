function [out,setupSWS]=FindSWS_XCorr_Each(inputXcorr,outplanes,setupdataprocessing,setupSWS,setupMakePlanes,dataDir)

% Estimate SWS based on input XCorr data
% inputXcorr contains the Xcorr data at each rotation angle, lateral position
% outplanes is a structure containing space-time planes for each acquistion rotation angle
% setupdataprocessing is a structure containing instructions for fitting
% setupSWS is a stucture containing information for fitting, and is modified and passed back out for record keeping
% setupMakePlane is a structure with plotting information
% dataDir is used to load B-mode information
% out is a structure contains SWS estimates at each rotation angle

% setup parameters
setupSWS=setTOFFitParams(setupdataprocessing,setupSWS,dataDir);

%% fit data for each angle
for iangle=1:length(outplanes)
    [out(iangle).SH]=fitSWStoXcorr(outplanes(iangle),inputXcorr.XCorrshifts_mspermm(iangle,:),inputXcorr.CorrCoeff(iangle,:),setupSWS.SH(iangle),setupSWS);
    out(iangle).SecondSW.speed=NaN;

end

%% Plot
if setupSWS.fignum
    set(0,'CurrentFigure',1);clf;

    nplot=2*setupdataprocessing.nangles;
    ncol=round((floor(sqrt(nplot*2.5)))/2)*2;
    nrow=ceil(nplot./ncol);
    for iangle=1:length(out)
        subplot(nrow,ncol,iangle*2-1);hold on
        imagesc(outplanes(iangle).tms,outplanes(iangle).latmm,outplanes(iangle).plane);
        caxis(setupMakePlanes.climsforplane)

        plot(out(iangle).SH.timesatlat,out(iangle).SH.latpos,'r')
        title([num2str(out(iangle).SH.speed,2) ' m/s'])
        tmp=gca;tmp.YDir='reverse';

        subplot(nrow,ncol,iangle*2);hold on;

        plot(outplanes(iangle).latmm(out(iangle).SH.ibelowcoeffthreshold),out(iangle).SH.Xcorrshifts(out(iangle).SH.ibelowcoeffthreshold),'.','Color',[.8 .5 .8]) % plot bad corrcoeff in pink
        plot(outplanes(iangle).latmm(~out(iangle).SH.ibelowcoeffthreshold),out(iangle).SH.Xcorrshifts(~out(iangle).SH.ibelowcoeffthreshold),'.','Color',[.7 .7 .7]) % plot good corrcoeff in grey
        plot(out(iangle).SH.latpos,out(iangle).SH.Xcorrshiftsfit,'.','Color',[.5 .8 .8]) % points that were sent to fit in blue
        plot(out(iangle).SH.latpos(out(iangle).SH.iTOFinlier),out(iangle).SH.Xcorrshiftsfit(out(iangle).SH.iTOFinlier),'k.') % plot points that are inliers for fit in black on top

        plot(out(iangle).SH.latpos,out(iangle).SH.shiftplot,'g-','LineWidth',1)
        title([num2str(setupdataprocessing.anglesDeg(iangle)) '^o'])
        ylim([-5 5])
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  setupSWS=setTOFFitParams(setupdataprocessing,setupSWS,dataDir)
if regexp(setupdataprocessing.SWSEstimationMethod,'.*FixedLat.*')
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_lat(\d+\.?\d?)to(\d+\.?\d?)outlierthres(\d+\.?\d*)Xcorrcoeffthres(\d+\.?\d*).*','tokens');
    latstarttmp=str2double(valstring{1}{1});
    latendtmp=str2double(valstring{1}{2});
    outlierthresval=str2double(valstring{1}{3});
    Xcorrcoeffthresholdval=str2double(valstring{1}{4});

    latstart(1:setupdataprocessing.nangles)=latstarttmp;
    latend(1:setupdataprocessing.nangles)=latendtmp;
    setupSWS.fittype='setlat';
elseif regexp(setupdataprocessing.SWSEstimationMethod,'.*AngleInformedLat.*')
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_Alonglat(\d+\.?\d?)to(\d+\.?\d?)Acrosslat(\d+\.?\d?)to(\d+\.?\d?)outlierthres(\d+\.?\d?)Xcorrcoeffthres(\d+\.?\d?)*','tokens');
    latstartAlong=str2double(valstring{1}{1});
    latendAlong=str2double(valstring{1}{2});
    latstartAcross=str2double(valstring{1}{3});
    latendAcross=str2double(valstring{1}{4});
    outlierthresval=str2double(valstring{1}{5});
    Xcorrcoeffthresholdval=str2double(valstring{1}{6});

    [latstart,latend]=setAngleInformedRange(latstartAlong,latendAlong,latstartAcross,latendAcross,setupdataprocessing,dataDir);
    setupSWS.fittype='setlat';
elseif regexp(setupdataprocessing.SWSEstimationMethod,'.*DynamicLat.*')
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_outlierthres(\d+\.?\d?)Xcorrcoeffthres(\d+\.?\d?).*minlatrange(\d+\.?\d?)Xcorrmeanerrorthres(\d+\.?\d?).*','tokens');
    outlierthresval=str2double(valstring{1}{1});
    Xcorrcoeffthresholdval=str2double(valstring{1}{2});
    minlatrange=str2double(valstring{1}{3});
    meanerrorthres=str2double(valstring{1}{4});

    latstart(1:setupdataprocessing.nangles)=0;
    latend(1:setupdataprocessing.nangles)=Inf;
    setupSWS.minlatrangelengthmm=minlatrange;
    setupSWS.meanerrorthres=meanerrorthres;
    setupSWS.fittype='dynamiclat';
end

for i=1:setupdataprocessing.nangles
    %set up SWSEstimationParams
    setupSWS.SH(i).minlatmm = latstart(i);
    setupSWS.SH(i).maxlatmm = latend(i);

    %other settings
    % copy global setting to where used
    setupSWS.SH(i).minspeedval = setupSWS.minspeed;
    setupSWS.SH(i).maxspeedval = setupSWS.maxspeed;
end

setupSWS.outlierthres=outlierthresval;
setupSWS.Xcorrcoeffthreshold=Xcorrcoeffthresholdval;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=fitSWStoXcorr(outplane,Xcorrshifts,Xcorrcoeffs, setupSWSindividual,setupSWS)

% make variable to manipulate
Xcorrshiftsfit=Xcorrshifts;

% remove xcorrs with bad correlation coefficents
ibelowthres=Xcorrcoeffs<setupSWS.Xcorrcoeffthreshold;
Xcorrshiftsfit(ibelowthres)=NaN;

% subselect lateral range
ilatposstart=find(outplane.latmm>=setupSWSindividual.minlatmm,1,'first');
ilatposend=find(outplane.latmm<=setupSWSindividual.maxlatmm,1,'last');
Xcorrshiftsfit(1:(ilatposstart-1))=NaN; % remove things outside of valid range
Xcorrshiftsfit((ilatposend+1):end)=NaN; % remove things outside of valid range
latpos=outplane.latmm; % keep record of NaN's with lateral position
latpos(1:(ilatposstart-1))=NaN; % remove things outside of valid range
latpos((ilatposend+1):end)=NaN; % remove things outside of valid range


if all(isnan(Xcorrshiftsfit))
    out.speed=NaN;
    out.latpos=latpos;
    out.shiftplot=NaN*ones(size(latpos));
    out.timesatlat=NaN*latpos;
    out.Xcorrshiftsfit=Xcorrshiftsfit;
    out.iTOFinlier=[];
    out.errorval=NaN;
    out.Xcorrshifts=Xcorrshifts;
    out.ibelowcoeffthreshold=ibelowthres;
    out.hittingedgeflag=0; %for plotting
    out.qualmetric=0;
    out.ttpqualmetric=NaN; %placeholder

    return
end

% find speed
switch setupSWS.fittype
    case 'setlat'
        out=findspeed(Xcorrshiftsfit,latpos,setupSWS,setupSWSindividual);
    case 'dynamiclat'
        out=dynamiclatrangefindspeed(Xcorrshiftsfit,latpos,setupSWS,setupSWSindividual);
end

%organize output
out.Xcorrshifts=Xcorrshifts;
out.ibelowcoeffthreshold=ibelowthres;
out.hittingedgeflag=0; %for plotting
out.qualmetric=0;
out.ttpqualmetric=out.errorval; %placeholder

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=findspeed(Xcorrshiftsfit,latpos,setupSWS,setupSWSindividual)

% setup fit: choose from set of shift options
Xcorrshiftsfit=reshape(Xcorrshiftsfit,[],1);
speedoptions=linspace(setupSWSindividual.minspeedval,setupSWSindividual.maxspeedval,500);
shiftoptions=1./speedoptions; %shifts in ms_mm
shiftoptions=reshape(shiftoptions,1,[]);

% calcuate difference to points
difftoval_squared=(shiftoptions-Xcorrshiftsfit).^2;

% max out outliers
difftoval_squared(difftoval_squared>setupSWS.outlierthres)=setupSWS.outlierthres;

%determine shift that leads to minimum error
summedsquarederror=sum(difftoval_squared,'omitnan');
[minval,imin]=min(summedsquarederror);
shiftvalue=shiftoptions(imin); % in delta(tms)/mm

%organize output
out.speed=1./shiftvalue;
out.latpos=latpos;
out.shiftplot=shiftvalue*ones(size(latpos));
out.timesatlat=(1/out.speed)*latpos;
out.Xcorrshiftsfit=Xcorrshiftsfit;
out.iTOFinlier=(Xcorrshiftsfit-shiftvalue).^2<setupSWS.outlierthres;
out.errorval=minval;

end
%%%%%%%%%%%%%%%%%%%%%%%%%

function out=dynamiclatrangefindspeed(Xcorrshiftsfit,latpos,setupSWS,setupSWSindividual)

% setup fit: choose from set of shift options
Xcorrshiftsfit=reshape(Xcorrshiftsfit,[],1);
speedoptions=linspace(setupSWSindividual.minspeedval,setupSWSindividual.maxspeedval,500);
shiftoptions=1./speedoptions; %shifts in ms_mm
shiftoptions=reshape(shiftoptions,1,[]);

%calculate difference to points
difftoval_squared=(shiftoptions-Xcorrshiftsfit).^2;

% max out outliers
difftoval_squared(difftoval_squared>setupSWS.outlierthres)=setupSWS.outlierthres;

% add dimension corresponding to start position for lateral range
% considered
difftoval_squared_lat=reshape(difftoval_squared,length(latpos),1,[]);
difftoval_squared_lat=repmat(difftoval_squared_lat,1,length(latpos),1);
[ilat,allilatstart]=meshgrid(1:length(latpos),1:length(latpos));
difftoval_squared_lat(repmat(allilatstart<ilat,[1,1,length(shiftoptions)]))=NaN; % set variable ilatstartpositions as second index

% cumsum to calculate for all combos of start, stop
summedsquarederror= cumsum(difftoval_squared_lat,1,'omitnan'); %sum over lateral position, now ilatend is first index

% find best shift value for each combo of start & stop options
[~,imin]=min(summedsquarederror,[],3,'linear');
[~,~,iminshift]=ind2sub(size(summedsquarederror),imin); % ilatend by ilatstart
shiftval=shiftoptions(iminshift); % ilatend by ilatstart, best shift value for given lat range

% now pick best lat range based on mean error from "best" shift value
% (already determined), only considering inliers
% expand Xcorrshifts so can compare to each "best" shift value for each lat
% range
Xcorrshiftsfittmp=repmat(Xcorrshiftsfit,1,length(latpos),length(latpos)); %ilat by ilatend by ilatstart
shiftval=repmat(reshape(shiftval,1,length(latpos),length(latpos)),length(latpos),1,1); %ilat by ilatend by ilatstart
[allilatend,ilat,allilatstart]=meshgrid(1:length(latpos),1:length(latpos),1:length(latpos));
shiftval(ilat<allilatstart)=NaN; % nan out data from before ilatstart
shiftval(ilat>allilatend)=NaN; % nan out data from after ilatend
iinlier=(shiftval-Xcorrshiftsfittmp).^2<setupSWS.outlierthres; % determine all inliers
Xcorrshiftsfittmp(~iinlier)=NaN; %nan out any outliers
sumsquarederror_inliersonly=squeeze(sum((shiftval-Xcorrshiftsfittmp).^2,1,'omitnan')); % recalculate sum error

% determine mean error
[allilatstart,allilatend]=meshgrid(1:length(latpos),1:length(latpos));
linelength=squeeze(allilatend-allilatstart);
meansquarederror_inliersonly=sumsquarederror_inliersonly./linelength;

% exclude anything that doesn't have the minimum lateral range
iminlatlength=ceil(setupSWS.minlatrangelengthmm./mean(diff(latpos)));
[allilatstart,allilatend]=meshgrid(1:length(latpos),1:length(latpos));
linelength=squeeze(allilatend-allilatstart);
meansquarederror_inliersonly(linelength<iminlatlength)=NaN;

% exclude anything that doesn't have the valid minimum number of points
minnumberinclude=floor(iminlatlength*.5); % at least 50% of points in minimum range have to be valid
numincluded=sum(~isnan(Xcorrshiftsfittmp),1);
meansquarederror_inliersonly(~(numincluded>minnumberinclude))=NaN;

% find 'best' lat range: longest 'great' error or, if no 'great' errors,
% lowest error
if ~sum(meansquarederror_inliersonly(:)<setupSWS.meanerrorthres) % if nothing below the threshold
    [~,ibestlatrange]=max(meansquarederror_inliersonly,[],[1 2],'linear'); % just use the max
else % otherwise, use longest line length that meets criteria
    allcheck=find(meansquarederror_inliersonly<setupSWS.meanerrorthres);
    [~,i3]=max(linelength(allcheck));
    ibestlatrange=allcheck(i3);
end

% refind best line with best lateral range
[ilatend,ilatstart]=ind2sub(size(meansquarederror_inliersonly),ibestlatrange);
summedsquarederror=sum(difftoval_squared(ilatstart:ilatend,:),'omitnan');
[minval,imin]=min(summedsquarederror);
finalshiftvalue=shiftoptions(imin); % in delta(tms)/mm

% update lateral range with points actually used
iinlier=find((Xcorrshiftsfit(ilatstart:ilatend)-finalshiftvalue).^2<setupSWS.outlierthres);
if ~isempty(iinlier)
    ilatstartupdated=iinlier(1)+ilatstart-1;ilatendupdated=iinlier(end)+ilatstart-1;
else
    ilatstartupdated=1;ilatendupdated=1;
end
% organize output
shiftplot=NaN(length(latpos),1);shiftplot(ilatstartupdated:ilatendupdated)=finalshiftvalue;
out.speed=1./finalshiftvalue;
out.latpos=latpos;
out.timesatlat=(1/out.speed)*latpos.*~isnan(shiftplot');
out.Xcorrshiftsfit=Xcorrshiftsfit;
out.shiftplot=shiftplot;
out.iTOFinlier=iinlier+ilatstart-1;
out.errorval=minval;

end
