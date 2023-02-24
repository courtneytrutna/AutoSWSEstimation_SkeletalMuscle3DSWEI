function [out,setupSWS]=FindSWS_TTP_Each(inputTOF,outplanes,setupdataprocessing,setupSWS,setupMakePlanes,dataDir)

% Estimate SWS at each rotation angle using the TTP method (all lateral range methods)
% inputTOF contains the TTP data at each rotation angle, lateral position
% outplanes is a structure containing space-time planes for each acquistion rotation angle
% setupdataprocessing is a structure containing instructions for fitting
% setupSWS is a stucture containing information for fitting, and is modified and passed back out for record keeping
% setupMakePlane is a structure with plotting information
% dataDir is used to load B-mode information
% out is a structure contains SWS estimates at each rotation angle

%% setup parameters
setupSWS=setTOFFitParams(setupdataprocessing,setupSWS,dataDir);

%% fit data for each angle
for iangle=1:length(outplanes)
    [out(iangle).SH]=fitSWStoTOF(outplanes(iangle),inputTOF(iangle,:),setupSWS.SH(iangle),setupSWS.outlierTOFthres,setupSWS.typeoffit);
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

        plot(out(iangle).SH.times,out(iangle).SH.latpos,'w')
        title([num2str(out(iangle).SH.speed,2) ' m/s'])
        tmp=gca;tmp.YDir='reverse';

        subplot(nrow,ncol,iangle*2);hold on;

        plot(outplanes(iangle).latmm,out(iangle).SH.TOF,'.','Color',[.7 .7 .7]) % plot all in grey
        plot(out(iangle).SH.latpos(out(iangle).SH.iTOFinlier),out(iangle).SH.TOFfit(out(iangle).SH.iTOFinlier),'k.') % plot points that used in black on top

        plot(out(iangle).SH.latpos,out(iangle).SH.times,'g-','LineWidth',1)
        title([num2str(setupdataprocessing.anglesDeg(iangle)) '^o, R^2:' num2str(out(iangle).SH.R2,2)])
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  setupSWS=setTOFFitParams(setupdataprocessing,setupSWS,dataDir)
if regexp(setupdataprocessing.SWSEstimationMethod,'.*FixedLat.*')
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_lat(\d+\.?\d?)to(\d+\.?\d?)originthres(\d+\.?\d?)outlierTOFthres(\d+\.?\d?).*','tokens');
    latstarttmp=str2double(valstring{1}{1});
    latendtmp=str2double(valstring{1}{2});
    originthresval=str2double(valstring{1}{3});
    outlierTOFthresval=str2double(valstring{1}{4});

    latstart(1:setupdataprocessing.nangles)=latstarttmp;
    latend(1:setupdataprocessing.nangles)=latendtmp;
    setupSWS.typeoffit='setlat';
elseif regexp(setupdataprocessing.SWSEstimationMethod,'.*AngleInformedLat.*')
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_Alonglat(\d+\.?\d?)to(\d+\.?\d?)Acrosslat(\d+\.?\d?)to(\d+\.?\d?)originthres(\d+\.?\d?)outlierTOFthres(\d+\.?\d?)*','tokens');
    latstartAlong=str2double(valstring{1}{1});
    latendAlong=str2double(valstring{1}{2});
    latstartAcross=str2double(valstring{1}{3});
    latendAcross=str2double(valstring{1}{4});
    originthresval=str2double(valstring{1}{5});
    outlierTOFthresval=str2double(valstring{1}{6});

    [latstart,latend]=setAngleInformedRange(latstartAlong,latendAlong,latstartAcross,latendAcross,setupdataprocessing,dataDir);
    setupSWS.typeoffit='setlat';
elseif regexp(setupdataprocessing.SWSEstimationMethod,'.*DynamicLat.*')
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_originthres(\d+\.?\d?)outlierTOFthres(\d+\.?\d?)minlatrange(\d+\.?\d?)R2thres(\d+\.?\d?).*','tokens');
    originthresval=str2double(valstring{1}{1});
    outlierTOFthresval=str2double(valstring{1}{2});
    minlatrangeval=str2double(valstring{1}{3});
    R2thresval=str2double(valstring{1}{4});

    latstart(1:setupdataprocessing.nangles)=0;
    latend(1:setupdataprocessing.nangles)=Inf;

    for i=1:setupdataprocessing.nangles
        setupSWS.SH(i).minlatrangelengthmm=minlatrangeval; %params specific to dyanmic range setting, so set here
        setupSWS.SH(i).R2thres=R2thresval;
    end
    setupSWS.typeoffit='dynamiclat';
end

for i=1:setupdataprocessing.nangles
    %set up SWSEstimationParams
    setupSWS.SH(i).minlatmm = latstart(i);
    setupSWS.SH(i).maxlatmm = latend(i);
    setupSWS.SH(i).originthresval=originthresval;

    %other settings
    setupSWS.SH(i).minspeedval = setupSWS.minspeed;
    setupSWS.SH(i).maxspeedval = setupSWS.maxspeed;
    setupSWS.outlierTOFthres=outlierTOFthresval;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=fitSWStoTOF(outplane,TOF,setupSWSindividual,outlierTOFthres,typeoffit)

ilatposstart=find(outplane.latmm>setupSWSindividual.minlatmm,1,'first');
ilatposend=find(outplane.latmm<setupSWSindividual.maxlatmm,1,'last');
TOFfit=TOF(ilatposstart:ilatposend); %only populate lateral-range valid regions
latpos=outplane.latmm(ilatposstart:ilatposend);

if all(isnan(TOFfit))
    out.speed=NaN;
    out.R2=NaN;
    out.latpos=latpos;
    out.times=nan(size(latpos));
    out.TOFfit=TOFfit;
    out.iTOFinlier=[];
    out.errorval=NaN;
    out.qualmetric= NaN;
    out.TOF=TOF;
    out.hittingedgeflag=0; %for plotting

    return
end

switch typeoffit
    case 'setlat'
        out=linearfit(TOFfit,latpos,setupSWSindividual,outlierTOFthres);
    case 'dynamiclat'
        out=dynamiclatlinearfit(TOFfit,latpos,setupSWSindividual,outlierTOFthres);

end
out.TOF=TOF;

out.hittingedgeflag=0; %for plotting

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=linearfit(TOFfit,latpos,setupSWS,outlierTOFthres)
% find best of predetermined options with outlier removal
% generate predetermined options
speedoptions=linspace(setupSWS.minspeedval,setupSWS.maxspeedval,500);
slopeoptions=1./speedoptions;slopeoptions=reshape(slopeoptions,1,[],1);
interceptoptions=linspace(-setupSWS.originthresval,setupSWS.originthresval,500);interceptoptions=reshape(interceptoptions,1,1,[]);
latpos=reshape(latpos,[],1);TOFfit=reshape(TOFfit,[],1);
lineoptions=latpos.*slopeoptions+interceptoptions;
% test all options
difftolines_squared=(lineoptions-TOFfit).^2;
% max out contriubtion of outliers
difftolines_squared(difftolines_squared>outlierTOFthres)=outlierTOFthres;
%difftolines_squared(isnan(difftolines_squared))=outlierTOFthres;
% pick best
summedsquarederror=sum(difftolines_squared,'omitnan');

[minerrorval,imin]=min(summedsquarederror(:));
[~,iminslope,iminintercept]=ind2sub(size(summedsquarederror),imin);

bestline(1)=slopeoptions(iminslope);
bestline(2)=interceptoptions(iminintercept);

%calc error vals for final, only considering inliers
lineval=lineoptions(:,iminslope,iminintercept);
iinlier=find((TOFfit-lineval).^2<outlierTOFthres);
SSres=sum((TOFfit(iinlier)-lineval(iinlier)).^2);
SStot=sum((TOFfit(iinlier)-mean(TOFfit(iinlier))).^2);

out.speed=1./bestline(1);
out.R2=1-(SSres/SStot);

out.latpos=latpos;
out.times=lineval;
out.TOFfit=TOFfit;
out.iTOFinlier=iinlier;
out.errorval=minerrorval;
out.qualmetric= out.errorval;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=dynamiclatlinearfit(TOFfit,latpos,setupSWS,outlierTOFthres)

speedoptions=linspace(setupSWS.minspeedval,setupSWS.maxspeedval,100);
interceptoptions=linspace(-setupSWS.originthresval,setupSWS.originthresval,25);

for i=1:3
    [out]=findbestlinefromoptions(TOFfit,latpos,setupSWS,outlierTOFthres,speedoptions,interceptoptions);
    speedoptions=ReduceSearchRangeAroundBest(speedoptions,out.speed,setupSWS.minspeedval,setupSWS.maxspeedval);
    interceptoptions=ReduceSearchRangeAroundBest(interceptoptions,out.intercept);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=findbestlinefromoptions(TOFfit,latpos,setupSWS,outlierTOFthres,speedoptions,interceptoptions)
% find best of predetermined options with outlier removal
% generate predetermined options

slopeoptions=1./speedoptions;slopeoptions=reshape(slopeoptions,1,[],1);
interceptoptions=reshape(interceptoptions,1,1,[]);
latpos=reshape(latpos,[],1);TOFfit=reshape(TOFfit,[],1);
lineoptions=latpos.*slopeoptions+interceptoptions;
% test all options
difftolines_squared=(lineoptions-TOFfit).^2;
% max out contriubtion of outliers
difftolines_squared(difftolines_squared>outlierTOFthres)=outlierTOFthres;

% determine minimum valid points
iminlatlength=ceil(setupSWS.minlatrangelengthmm./mean(diff(latpos)));
minnumberinclude=floor(iminlatlength*.5); % at least 50% of points in minimum range have to be valid

difftolines_squared_lat=reshape(difftolines_squared,length(latpos),1,length(slopeoptions),length(interceptoptions));
difftolines_squared_lat=repmat(difftolines_squared_lat,1,length(latpos),1,1);
[ilat,allilatstart]=meshgrid(1:length(latpos),1:length(latpos));
difftolines_squared_lat(repmat(allilatstart<ilat,1,1,length(slopeoptions),length(interceptoptions)))=NaN; % set variable ilatstartpositions as second index

summedsquarederror=cumsum(difftolines_squared_lat,1,'omitnan'); % sum over lateral position, now ilatend is first index
[~,imin]=min(summedsquarederror,[],[3 4],'linear'); % find best line for each combo of start & stop options

[~,~,iminslope,iminintercept]=ind2sub(size(summedsquarederror),imin);

bestslope=slopeoptions(iminslope);bestslope=reshape(bestslope,1,length(latpos),length(latpos));
bestintercept=interceptoptions(iminintercept);bestintercept=reshape(bestintercept,1,length(latpos),length(latpos));
lineval=latpos.*bestslope+bestintercept; %lat by ilatend by ilatstart
[allilatend,ilat,allilatstart]=meshgrid(1:length(latpos),1:length(latpos),1:length(latpos));
lineval(ilat<allilatstart)=NaN; % nan out data from before ilatstart
lineval(ilat>allilatend)=NaN; % nan out data from after ilatend

TOFfittmp=repmat(TOFfit,1,length(latpos),length(latpos));

iinlier=(TOFfittmp-lineval).^2<outlierTOFthres; % determine all inliers
TOFfittmp(~iinlier)=NaN;
SSres=sum(((TOFfittmp-lineval)).^2,1,'omitnan');
SStot=sum((TOFfittmp-mean(TOFfittmp,1,'omitnan')).^2,1,'omitnan');
R2=squeeze(1-(SSres./SStot)); % ilatend by ilatstart

numincluded=sum(~isnan(TOFfittmp),1);
R2(~(numincluded>minnumberinclude))=NaN;

[allilatstart,allilatend]=meshgrid(1:length(latpos),1:length(latpos));
linelength=squeeze(allilatend-allilatstart);
R2(linelength<iminlatlength)=NaN;

if ~sum(R2(:)>setupSWS.R2thres) % if nothing above the threshold
    [~,ibestlatrange]=max(R2,[],[1 2],'linear'); % just use the max
else % otherwise, use longest line length that meets criteria
    allcheck=find(R2>setupSWS.R2thres);
    [~,i3]=max(linelength(allcheck));
    ibestlatrange=allcheck(i3);
end

% refind best line with best lateral range
[ilatend,ilatstart]=ind2sub(size(R2),ibestlatrange);
summedsquarederror=sum(difftolines_squared(ilatstart:ilatend,:,:),'omitnan');
[minerrorval,imin]=min(summedsquarederror,[],[1 2 3],'linear');
[~,iminslope,iminintercept]=ind2sub(size(summedsquarederror),imin);
bestline(1)=slopeoptions(iminslope);
bestline(2)=interceptoptions(iminintercept);
%calc error vals for final, only considering inliers
lineval=lineoptions(:,iminslope,iminintercept);
TOFfittmp=TOFfit(ilatstart:ilatend); linevaltmp=lineval(ilatstart:ilatend);
iinlier=find((TOFfittmp-linevaltmp).^2<outlierTOFthres);
SSres=sum((TOFfittmp(iinlier)-linevaltmp(iinlier)).^2);
SStot=sum((TOFfittmp(iinlier)-mean(TOFfittmp(iinlier))).^2);
% adjust what ilatstart and ilatend were based on which points actually got
% included
ilatstartupdated=iinlier(1)+ilatstart-1;ilatendupdated=iinlier(end)+ilatstart-1;

out.speed=1./bestline(1);
out.intercept=bestline(2);
out.R2=1-(SSres/SStot);

lineval=zeros(size(lineval)); lineval(ilatstartupdated:ilatendupdated)=linevaltmp(iinlier(1):iinlier(end));

out.latpos=latpos;
out.times=lineval;
out.TOFfit=TOFfit;
out.iTOFinlier=iinlier+ilatstart-1;
out.errorval=minerrorval;
out.qualmetric= out.errorval;
end