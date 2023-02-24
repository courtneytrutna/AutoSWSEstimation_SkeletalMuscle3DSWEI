function [out3DSWS,setupSWS]=FindSWS_TTP_SetLatRange_3D(inputTTP,outplanes,setupdataprocessing,setupSWS,setupMakePlanes,dataDir)

% Estimate SWS based on input TTP data, fitting entire 3D dataset at once
% inputTOF contains the TTP data at each rotation angle, lateral position
% outplanes is a structure containing space-time planes for each acquistion rotation angle
% setupdataprocessing is a structure containing instructions for fitting
% setupSWS is a stucture containing information for fitting, and is modified and passed back out for record keeping
% setupMakePlane is a structure with plotting information
% dataDir is used to load B-mode information
% out3DSWS is a structure contains SWS estimates for the 3D-SWEI acquisition

% set method for surface fitting
setupSWS.fit3Dmethod='direct';
% ransac method is also coded, but direct fit is method used in paper and
% recommended

%% setup parameters
setupSWS=setTOFFitParams(setupdataprocessing,setupSWS,dataDir);

%% Fit all TOF simultaneously
[out]=fitSWStoTOFsurface(outplanes(1),inputTTP,setupSWS,setupdataprocessing.anglesDeg);

%% Oraganize Output
% Final outvalues: consistent format across all methods (cPar,cPerp) plus
% holding on to information of interest

out3DSWS.cPar=out.cPar;
out3DSWS.cPerp=out.cPerp;
out3DSWS.phiRot=out.phirot;

out=rmfield(out,{'cPar','cPerp','phirot'});

out3DSWS.EachSWSEstimate=out;

%% Plotting
if setupSWS.fignum

    figure(2);clf;
    %plot3DTTP(expanglesinput,latmm,squeeze(tmssurface(:,:,icparmin,icperpmin,iphirotmin)-mediandiff(:,:,icparmin,icperpmin,iphirotmin)),'surf')
    plot3Dhelper(setupdataprocessing.anglesDeg,outplanes(1).latmm,out3DSWS.EachSWSEstimate.surftms,'surf');

    hold on;
    tmpplotinlier=out3DSWS.EachSWSEstimate.TOFtofit;tmpplotinlier(~out3DSWS.EachSWSEstimate.iTOFinlier)=NaN;
    plot3Dhelper(setupdataprocessing.anglesDeg,outplanes(1).latmm,tmpplotinlier,'k.');
    tmpplotoutlier=out3DSWS.EachSWSEstimate.TOFtofit;tmpplotoutlier(out3DSWS.EachSWSEstimate.iTOFinlier)=NaN;
    plot3Dhelper(setupdataprocessing.anglesDeg,outplanes(1).latmm,tmpplotoutlier,'g.');

    plot3Dhelper([],[],[],'grid');
    plot3(cosd(out3DSWS.phiRot)*[0 20],sind(out3DSWS.phiRot)*[0 20],[0 0],'b-','LineWidth',2)

    colorline=[.6 .6 .6];
    text(20,0,0,'0^o','Color',colorline)

    for iz=[5 10 15]
        plot3([-20 20],[20 20],[iz iz],'-','Color',colorline)
        plot3([-20 -20],[-20 20],[iz iz],'-','Color',colorline)
    end

    for ir=[4 8 12 16 20]
        plot3(ir*cos(linspace(0,2*pi,100)),ir*sin(linspace(0,2*pi,100)),0*linspace(0,2*pi,100),'-','Color',colorline)
        text(1,-ir+1,0,[num2str(ir) ' mm'],'Color',colorline)
    end

    save3Dplot=1;
    if save3Dplot
        [savefolder,SWSsettingsname]=GenerateSaveFileName(dataDir,setupdataprocessing);
        filenamegif=[savefolder 'TOF3DGif' SWSsettingsname '.gif'];
        tmp=gca;tmp.View= [70 18];% [104 24]

        tmp.XTick=[];
        tmp.YTick=[];
        tmp.ZTick=[0 5 10 15];

        V1=linspace(70,125,50);
        V2=linspace(18,18,50);
        figure(2);
        f1=gcf;f1.Color=[1 1 1];f1.Position=[575 875 1200 900];

        for i=[1:50 50*ones([1,50])]
            tmp.View=[V1(i) V2(i)];
            gifify(f1,i,filenamegif)
        end
        for i=[50:-1:1 ones([1,10])]
            tmp.View=[V1(i) V2(i)];
            gifify(f1,0,filenamegif)
        end
    end


    % Plot
    set(0,'CurrentFigure',1);clf;

    nplot=2*setupdataprocessing.nangles;
    ncol=round((floor(sqrt(nplot*2.5)))/2)*2;
    nrow=ceil(nplot./ncol);
    ymax=max(inputTTP(:));
    for iangle=1:length(outplanes)
        subplot(nrow,ncol,iangle*2-1);hold on
        imagesc(outplanes(iangle).tms,outplanes(iangle).latmm,outplanes(iangle).plane);
        caxis(setupMakePlanes.climsforplane)

        plot(out3DSWS.EachSWSEstimate.surftms(iangle,:),outplanes(iangle).latmm,'w')
        title([num2str(out3DSWS.EachSWSEstimate.cvalsall(iangle),2) ' m/s'])
        tmp=gca;tmp.YDir='reverse';

        subplot(nrow,ncol,iangle*2);hold on;

        plot(outplanes(iangle).latmm,inputTTP(iangle,:),'.','Color',[.7 .7 .7]) % plot all in grey
        plot(outplanes(iangle).latmm,tmpplotinlier(iangle,:),'k.')

        plot(outplanes(iangle).latmm,out3DSWS.EachSWSEstimate.surftms(iangle,:),'g-','LineWidth',1)
        xlim([outplanes(iangle).latmm(1) outplanes(iangle).latmm(end)]);ylim([0 ymax])
    end

    % move everything over a bit
    tmp=figure(1);
    tmp=findobj(tmp,'Type','Axes');
    for i=1:length(tmp);tmp(i).Position(1)=tmp(i).Position(1)-.1;end

    % add polar plot in top left
    polax=polaraxes;
    polax.Position=[.78 .65 .25 .25];
    hold on
    rlim([0 8])
    polarplot(setupdataprocessing.anglesDeg*pi/180,out3DSWS.EachSWSEstimate.cvalsall,'k-') %plot SWS data from fit

    polarplot([1 1]*out3DSWS.phiRot.*pi/180,[0 6],'k--')

    title({['ellipse cost function val: ' num2str(out3DSWS.EachSWSEstimate.errorval,3) ]})

    %plot values
    ax_vals=axes;
    ax_vals.Position=[.85 .4 .1 .15];
    offset1=0;
    offset2=0.1;
    plot((1:2)-offset1-offset2,[out3DSWS.cPar,out3DSWS.cPerp],'k*')
    legend('SH fit','Location','southoutside')
    tmp=gca;tmp.XTick=[1,2];tmp.XTickLabel={'c_{Par}','c_{Perp}'}; %label axis
    xlim([.5 2.5]);ylim([0 10])
    ylabel('m/s')
    title({['SWS_L= ' num2str(out3DSWS.cPar,3) '; SWS_T= ' num2str(out3DSWS.cPerp,3)]})
    ax_vals.Position=[.85 .4 .1 .15]; %reset after legend

    % plot angles
    ax_deg=axes;
    ax_deg.Position=[.85 .1 .1 .15];
    plot(1-offset1-offset2,out3DSWS.phiRot,'k*')
    hold on
    if isfield(setupSWS,'bmoderotangle')
        plot(1-offset1,setupSWS.bmoderotangle,'r*')
    end

    hold off
    tmp=gca;tmp.XTick=[];
    xlim([.5 1.5]);ylim([-90 90])
    ylabel('Rot Angle (degrees)')
    legendstrings={'SH fit'};
    if isfield(setupSWS,'bmoderotangle')
        legendstrings{length(legendstrings)+1}='Bmode';
    end
    legend(legendstrings,'Location','southoutside')
    title(['Rot =  ' num2str(out3DSWS.phiRot) '^o, Tilt = 0^o'])
    ax_deg.Position=[.85 .1 .1 .15]; %reset after legend

    tmp=regexp(dataDir,'/','split');titlestring=tmp{end};titlestring=replace(titlestring,'_',' ');

    sgtitle(titlestring)

    [savefolder,SWSsettingsname]=GenerateSaveFileName(dataDir,setupdataprocessing);
    print('-dpng',[savefolder 'SpaceTimeAndEllipsePlots' SWSsettingsname '.png'])

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

elseif regexp(setupdataprocessing.SWSEstimationMethod,'.*AngleInformedLat.*')
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_Alonglat(\d+\.?\d?)to(\d+\.?\d?)Acrosslat(\d+\.?\d?)to(\d+\.?\d?)originthres(\d+\.?\d?)outlierTOFthres(\d+\.?\d?)*','tokens');
    latstartAlong=str2double(valstring{1}{1});
    latendAlong=str2double(valstring{1}{2});
    latstartAcross=str2double(valstring{1}{3});
    latendAcross=str2double(valstring{1}{4});
    originthresval=str2double(valstring{1}{5});
    outlierTOFthresval=str2double(valstring{1}{6});


    [latstart,latend,bmoderotangle]=setAngleInformedRange(latstartAlong,latendAlong,latstartAcross,latendAcross,setupdataprocessing,dataDir);
    setupSWS.bmoderotangle=bmoderotangle;
end

for i=1:setupdataprocessing.nangles
    %set up SWSEstimationParams
    setupSWS.SH(i).minlatmm = latstart(i);
    setupSWS.SH(i).maxlatmm = latend(i);

    %other settings
end
setupSWS.outlierTOFthres=outlierTOFthresval;
setupSWS.maxspeedval=setupSWS.maxspeed;
setupSWS.minspeedval=setupSWS.minspeed;
setupSWS.origintmsthreshold=originthresval;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=fitSWStoTOFsurface(outplanes,TOF,setupSWS,expangles)

TOFfit=nan(size(TOF)); %initialize as nans
latmap=zeros(length(expangles),length(outplanes.latmm));
for iang=1:length(expangles)
    ilatposstart=find(outplanes.latmm>setupSWS.SH(iang).minlatmm,1,'first');
    ilatposend=find(outplanes.latmm<setupSWS.SH(iang).maxlatmm,1,'last');
    TOFfit(iang,ilatposstart:ilatposend)=TOF(iang,ilatposstart:ilatposend); %only populate lateral-range valid regions
    latmap(iang,ilatposstart:ilatposend)=1;
end

if strcmp(setupSWS.fit3Dmethod,'ransac')
    out=ransacsurfacefit(TOFfit,setupSWS,outplanes,expangles,latmap);
elseif strcmp(setupSWS.fit3Dmethod,'direct')
    out=fitTOFsurface(TOFfit,expangles,outplanes.latmm,latmap,setupSWS);
end

out.hittingedgeflag=0; %for plotting
out.tofqualmetric=0;
out.latpos=outplanes.latmm;
out.ilatinclude=latmap;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=ransacsurfacefit(TOFfit,setupSWS,outplanes,expangles,latmap)

niter=100;
besterrorachieved=length(latmap(:))*100; % this is a terrible error, and will almost certainly be passed.
currentbestsurface=[0 0];

rng('default'); % reset for reproducibility
%warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale') %supress commmon warnings
%warning('off','MATLAB:polyfit:PolyNotUnique')
for iiterations=1:niter
    itest=randi(length(TOFfit(:)),[50,1]); % chose X number of sample points to fit
    [itest1,itest2]=ind2sub(size(TOFfit),itest);
    TOFfittmp=NaN(size(TOFfit));
    for i=1:length(itest1)
        TOFfittmp(itest1(i),itest2(i))=TOFfit(itest1(i),itest2(i)); % include only those sample points
    end

    % fit SWS surface
    [tmpsurface]=fitTOFsurface(TOFfittmp,expangles,outplanes.latmm,latmap,setupSWS);
    % limits on surface fit are inside of fitTOF surface

    % only include inliers:
    iinlier=find((TOFfit-tmpsurface.surftms).^2<setupSWS.outlierTOFthres);

    %refit to all inlier points
    [iinlier1,iinlier2]=ind2sub(size(TOFfit),iinlier);
    TOFfittmp=NaN(size(TOFfit));
    for i=1:length(itest1)
        TOFfittmp(iinlier1(i),iinlier2(i))=TOFfit(iinlier1(i),iinlier2(i)); % include only those sample points
    end
    [tmpsurface]=fitTOFsurface(TOFfittmp,expangles,outplanes.latmm,latmap,setupSWS);

    %RANSAC quality metric
    currenterrormetric=tmpsurface.errorval; % for a given acq, nan number should be constant, so don't need to worry about reducing value by not including nans

    if currenterrormetric<besterrorachieved
        besterrorachieved=currenterrormetric;
        currentbestsurface=tmpsurface;
    end

end

out=currentbestsurface;

end


%%%%%%%%%%%%%%%%%%%%%%%%%
function outsurface=fitTOFsurface(TOFfit,expanglesinput,latmm,latmap,setupSWS)

% order of dimensions:
% angles x lateral x cPar x cPerp x phiRot

cParvals  = (setupSWS.minspeed*2):0.5:setupSWS.maxspeed;      % values to try "parallel" to fibers
cPerpvals = (setupSWS.minspeed):0.5:(setupSWS.maxspeed*.5);   % values to try "perpendicular" to fibers
phirotvals= -90:10:80;      % rotation angles to try
tshiftParvals=-setupSWS.origintmsthreshold:.5:setupSWS.origintmsthreshold;      % how many much time (in ms) to give as flux
tshiftPerpvals=-setupSWS.origintmsthreshold:.5:setupSWS.origintmsthreshold;      % how many much time (in ms) to give as flux


for i=1:3

    [outsurface]=calcBestSurfaceFromOptions(TOFfit,expanglesinput,latmm,latmap,setupSWS,cParvals,cPerpvals,phirotvals,tshiftParvals,tshiftPerpvals);

    cParvals=ReduceSearchRangeAroundBest(cParvals,outsurface.cPar);
    cPerpvals=ReduceSearchRangeAroundBest(cPerpvals,outsurface.cPerp);
    phirotvals=ReduceSearchRangeAroundBest(phirotvals,outsurface.phirot);
    tshiftParvals=ReduceSearchRangeAroundBest(tshiftParvals,outsurface.tshiftPar,-setupSWS.origintmsthreshold,setupSWS.origintmsthreshold);
    tshiftPerpvals=ReduceSearchRangeAroundBest(tshiftPerpvals,outsurface.tshiftPerp,-setupSWS.origintmsthreshold,setupSWS.origintmsthreshold);
end


end
%%%%%%%%%%%
function [outsurface]=calcBestSurfaceFromOptions(TOFfit,expanglesinput,latmm,latmap,setupSWS,cParvals,cPerpvals,phirotvals,tshiftParvals,tshiftPerpvals)

expanglesreshape(:,1,1,1,1,1)=single(expanglesinput);
latmmreshape(1,:,1,1,1,1)=single(latmm);
cParreshape(1,1,:,1,1,1)=single(cParvals);
cPerpreshape(1,1,1,:,1,1)=single(cPerpvals);
phirotreshape(1,1,1,1,:,1)=single(phirotvals);

% V=sqrt((muL * muT) / (muL*sin^2(psi)+muT*cos^2(psi))/rho)
% muL=cL^2*rho; muT=cT^2*rho
% psi=experimentalangle-ellipserotationangle
% V=sqrt((cL^2*rho * cT^2*rho) / (cL^2*rho*sin^2(psi)+cT^2*rho*cos^2(psi))/rho)
% rho's cancel out
% V=sqrt((cL^2 * cT^2) / (cL^2*sin^2(psi)+cT^2*cos^2(psi)))

cvalsall=sqrt((cParreshape.^2.*cPerpreshape.^2)./(cParreshape.^2.*(sind(expanglesreshape-phirotreshape).^2)+cPerpreshape.^2.*cosd(expanglesreshape-phirotreshape).^2));
clear cParreshape cParreshape
cvalsall=repmat(cvalsall,[1,size(TOFfit,2),1,1,1,1]); % all combos of speeds as a function of angle)

latmmmatrix=repmat(latmmreshape,[size(cvalsall,1),1,size(cvalsall,3),size(cvalsall,4),size(cvalsall,5),size(cvalsall,6)]);

tofsurface=latmmmatrix./cvalsall; %time of wavefront at every lat position, angle, for all combos of speeds
clear latmmmatrix latmmreshape

% calculate difference of points to all possible surfaces
TOFfitall=repmat(single(TOFfit),[1,1,size(tofsurface,3),size(tofsurface,4),size(tofsurface,5)]);
difftosurface=tofsurface-TOFfitall;
clear TOFfitall

% we're okay with 'wave didn't start propogating yet' shifts, so
% allow for set of shifts as additional fit param
tmptocorrectneg=10;
tshiftParreshape(1,1,1,1,1,:,1)=single(tshiftParvals)+tmptocorrectneg; %add 10 to keep track of negative, and then subtract off later
tshiftPerpreshape(1,1,1,1,1,1,:)=single(tshiftPerpvals)+tmptocorrectneg;

tshiftall=sqrt((tshiftParreshape.^2.*tshiftPerpreshape.^2)./(tshiftParreshape.^2.*(sind(expanglesreshape-phirotreshape).^2)+tshiftPerpreshape.^2.*(cosd(expanglesreshape-phirotreshape)).^2));
tshiftall=tshiftall-tmptocorrectneg; % this creates ellipse like shapes with positive and negative range

difftosurface_shift_squared=(difftosurface+tshiftall).^2;
clear expanglesreshape phirotreshape


% max out the contribution of outliers for each surface if doing a direct
% fit; if doing a RANSAC fit skip this step
if strcmp(setupSWS.fit3Dmethod,'direct')
    difftosurface_shift_squared(difftosurface_shift_squared>setupSWS.outlierTOFthres)=setupSWS.outlierTOFthres;
end

% the sum up the squared remaining differences
summeddiff=sum(difftosurface_shift_squared,[1,2],'omitnan')./sum(~isnan(difftosurface_shift_squared),[1 2]);

% and find minimum
[errorval,imin]=min(summeddiff(:));
[~,~,icparmin,icperpmin,iphirotmin,itshiftminpar,itshiftminperp]=ind2sub(size(summeddiff),imin);

% output relevant values
tmssurfaceplot=squeeze(tofsurface(:,:,icparmin,icperpmin,iphirotmin))+squeeze(tshiftall(:,:,:,:,iphirotmin,itshiftminpar,itshiftminperp)); %choose best surface
tmssurfaceplot(~latmap)=NaN; % only plot on considered lateral ranges

outsurface.cPar=cParvals(icparmin);
outsurface.cPerp=cPerpvals(icperpmin);
outsurface.phirot=phirotvals(iphirotmin);
outsurface.tshiftPar=tshiftParvals(itshiftminpar);
outsurface.tshiftPerp=tshiftPerpvals(itshiftminperp);
outsurface.surftms=tmssurfaceplot;
outsurface.cvalsall=squeeze(cvalsall(:,1,icparmin,icperpmin,iphirotmin)); %only need at one lateral position,
outsurface.TOFtofit=TOFfit;
outsurface.iTOFinlier=(TOFfit-tmssurfaceplot).^2<setupSWS.outlierTOFthres;
outsurface.errorval=errorval;

if outsurface.cPar<outsurface.cPerp % if along>across, flip values
    if outsurface.phirot>=0
        outsurface.phirot=outsurface.phirot-90;
    else
        outsurface.phirot=outsurface.phirot+90;
    end
    tmphold=outsurface.cPar;
    outsurface.cPar=outsurface.cPerp;
    outsurface.cPerp=tmphold;
end
end

%%%%%
