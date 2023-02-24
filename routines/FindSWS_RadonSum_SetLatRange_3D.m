function [out3DSWS,setupSWS]=FindSWS_RadonSum_SetLatRange_3D(outplanes,setupdataprocessing,setupSWS,setupMakePlanes,dataDir)

% estimate SWS with radon sum method over entire 3D dataset at once
%
% outplanes is a structure containing space-time planes for each acquistion rotation angle
% setupdataprocessing is a structure containing instructions for fitting
% setupSWS is a stucture containing information for fitting, and is modified and passed back out for record keeping
% setupMakePlane is a structure with plotting information
% dataDir is used to load B-mode information and for labeling and saving plots
% out is a structure contains SWS estimates at each rotation angle


%% setup parameters
setupSWS=setRadonSumParams(setupdataprocessing,setupSWS,dataDir);

%% make 3D dataset
[data3Dfit,latmap,data3Dfull]=make3DDatasetForRadonSum(outplanes,setupSWS);

%% Fit 3D radon sum
[out]=fit3DRadonSum(data3Dfit,latmap,setupSWS,setupdataprocessing.anglesDeg,outplanes(1).latmm,outplanes(1).tms);

%% Oraganize Output
% Final outvalues: consistent format across all methods (cPar,cPerp) plus
% holding on to information of interest

out3DSWS.cPar=out.cPar;
out3DSWS.cPerp=out.cPerp;
out3DSWS.phiRot=out.phirot;

out=rmfield(out,{'cPar','cPerp','phirot'});

out3DSWS.EachSWSEstimate=out;

%% Plot Output
if setupSWS.fignum
    Plot3DRadonSum(out3DSWS,setupdataprocessing,data3Dfit,data3Dfull,outplanes,setupSWS,setupMakePlanes,dataDir)
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function  setupSWS=setRadonSumParams(setupdataprocessing,setupSWS,dataDir)

if regexp(setupdataprocessing.SWSEstimationMethod,'.*FixedLat.*')
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_lat(\d+\.?\d?)to(\d+\.?\d?)originthres(\d+\.?\d?).*','tokens');
    latstarttmp=str2double(valstring{1}{1});
    latendtmp=str2double(valstring{1}{2});
    originthresval=str2double(valstring{1}{3});

    latstart(1:setupdataprocessing.nangles)=latstarttmp;
    latend(1:setupdataprocessing.nangles)=latendtmp;

elseif regexp(setupdataprocessing.SWSEstimationMethod,'.*AngleInformedLat.*')

    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_Alonglat(\d+\.?\d?)to(\d+\.?\d?)Acrosslat(\d+\.?\d?)to(\d+\.?\d?)originthres(\d+\.?\d?)*','tokens');
    latstartAlong=str2double(valstring{1}{1});
    latendAlong=str2double(valstring{1}{2});
    latstartAcross=str2double(valstring{1}{3});
    latendAcross=str2double(valstring{1}{4});
    originthresval=str2double(valstring{1}{5});

    [latstart,latend,bmoderotangle]=setAngleInformedRange(latstartAlong,latendAlong,latstartAcross,latendAcross,setupdataprocessing,dataDir);
    setupSWS.bmoderotangle=bmoderotangle;

else
    error('unimplemented method to set Radon Sum SWS Params')
end


for i=1:setupdataprocessing.nangles
    setupSWS.SH(i).minlatmm = latstart(i);
    setupSWS.SH(i).maxlatmm = latend(i);

    %setupSWS.SH(i).maxspeed = setupSWS.maxspeed;
    %setupSWS.SH(i).minspeed = setupSWS.minspeed;

    % leave these incase want later, but make super general because somewhat
    % redundant with above
    setupSWS.SH(i).tmsminstarttime=0; % require wave doesn't start in negative time
    setupSWS.SH(i).tmsmaxstart=100;
    setupSWS.SH(i).tmsmaxend=100;

end
setupSWS.origintmsthreshold=originthresval;
setupSWS.maxspeed = setupSWS.maxspeed;
setupSWS.minspeed = setupSWS.minspeed;


end

%%%%%%%%%%%%%%%%%%%%
function [dataset3Dtrimmed,latmap,dataset3Dfull]=make3DDatasetForRadonSum(outplanes,setupSWS)
dataset3Dtrimmed=nan(length(outplanes),size(outplanes(1).plane,1),size(outplanes(1).plane,2)); %initialize as nans
dataset3Dfull=nan(length(outplanes),size(outplanes(1).plane,1),size(outplanes(1).plane,2)); %initialize as nans
latmap=zeros(length(outplanes),length(outplanes(1).latmm));
for iang=1:length(outplanes)
    ilatposstart=find(outplanes(iang).latmm>setupSWS.SH(iang).minlatmm,1,'first');
    ilatposend=find(outplanes(iang).latmm<setupSWS.SH(iang).maxlatmm,1,'last');
    dataset3Dtrimmed(iang,ilatposstart:ilatposend,:)=outplanes(iang).plane(ilatposstart:ilatposend,:); %only populate lateral-range valid regions
    latmap(iang,ilatposstart:ilatposend)=1;
    dataset3Dfull(iang,:,:)=outplanes(iang).plane;
end

end

%%%%%%%%%%%%%%%%%%%%%

function [out]=fit3DRadonSum(data3D,latmap,setupSWS,expanglesinput,latmm,tms)

normalizationflag=1;
if normalizationflag==1
    medianperangle=median(data3D(:,:,:),[2 3],'omitnan');
    data3D=data3D-medianperangle;
elseif normalizationflag==2
    maxvals=max(data3D,[],[2 3]);
    data3D=data3D./maxvals;
end

cParvals  = (setupSWS.minspeed*2):0.5:setupSWS.maxspeed;      % values to try "parallel" to fibers
cPerpvals = (setupSWS.minspeed):0.5:(setupSWS.maxspeed*.5);   % values to try "perpendicular" to fibers
phirotvals= -90:10:80;      % rotation angles to try
tshiftvals=-setupSWS.origintmsthreshold:.5:setupSWS.origintmsthreshold;      % how many much time (in ms) to give as flux

for i=1:3
    [out]=FindBestRadonSumEllipseFromOptions(data3D,latmap,expanglesinput,latmm,tms,cParvals,cPerpvals,phirotvals,tshiftvals);

    cParvals=ReduceSearchRangeAroundBest(cParvals,out.cPar,setupSWS.minspeed,setupSWS.maxspeed);
    cPerpvals=ReduceSearchRangeAroundBest(cPerpvals,out.cPerp,setupSWS.minspeed,setupSWS.maxspeed);
    phirotvals=ReduceSearchRangeAroundBest(phirotvals,out.phirot);
    tshiftvals=ReduceSearchRangeAroundBest(tshiftvals,out.tshift,-setupSWS.origintmsthreshold,setupSWS.origintmsthreshold);
end

end

function [outRadonEllipse]=FindBestRadonSumEllipseFromOptions(data3D,latmap,expanglesinput,latmm,tms,cParvals,cPerpvals,phirotvals,tshiftvals)

expanglesreshape(:,1,1,1,1,1)=single(expanglesinput);
latmmreshape(1,:,1,1,1,1)=single(latmm);
cParreshape(1,1,:,1,1,1)=single(cParvals);
cPerpreshape(1,1,1,:,1,1)=single(cPerpvals);
phirotreshape(1,1,1,1,:,1)=single(phirotvals);
tshiftvalsreshape(1,1,1,1,1,:)=single(tshiftvals);

cvalsall=sqrt((cParreshape.^2.*cPerpreshape.^2)./(cParreshape.^2.*(sind(expanglesreshape-phirotreshape).^2)+cPerpreshape.^2.*cosd(expanglesreshape-phirotreshape).^2));
cvalsall=repmat(cvalsall,[1,size(data3D,2),1,1,1,1]); % all combos of speeds as a function of angle)

latmmmatrix=repmat(latmmreshape,[size(cvalsall,1),1,size(cvalsall,3),size(cvalsall,4),size(cvalsall,5),size(cvalsall,6)]);

twavefront=latmmmatrix./cvalsall+tshiftvalsreshape; %time of wavefront at every lat position, angle, for all combos of speeds
it=(twavefront-tms(1))./median(diff(tms)); %which ideal temporal index this wavefront time corresponds to
it_pre=floor(it); %time step before ideal
interpFrac=1-(it-it_pre); % interp fraction
% used to be like: sum(plane(it_pre).*fracs+plane(it_pre+1).*(1-fracs))
for icpar=1:size(it_pre,3)
    for icerp=1:size(it_pre,4)
        for iphirot=1:size(it_pre,5)
            for itshift=1:size(it_pre,6)
                it_pre_onemat=squeeze(it_pre(:,:,icpar,icerp,iphirot,itshift));
                interpFrac_onemat=squeeze(interpFrac(:,:,icpar,icerp,iphirot,itshift));
                ilatindices=ones(size(data3D,1),1)*(1:length(latmm));
                ianglesindices=(1:size(data3D,1))'*ones(1,size(data3D,2));

                % some speeds are too fast! Set to 1 so runs but then make interpfactor nan so not include
                interpFrac_onemat(it_pre_onemat<1)=NaN;
                it_pre_onemat(it_pre_onemat<1)=1;
                % ditto for speeds that are too slow
                interpFrac_onemat(it_pre_onemat>size(data3D,3))=NaN;
                it_pre_onemat(it_pre_onemat>(size(data3D,3)-1))=size(data3D,3)-1; %need -1 because of +1 for interp

                normalizationfactor=sum(~isnan(interpFrac_onemat(:))); %normalize for dropping points because too fast or slow

                indices1test = sub2ind(size(data3D),ianglesindices(:),ilatindices(:),it_pre_onemat(:));
                indices2test = sub2ind(size(data3D),ianglesindices(:),ilatindices(:),it_pre_onemat(:)+1);

                radonout(icpar,icerp,iphirot,itshift)=sum(data3D(indices1test).*interpFrac_onemat(:)+data3D(indices2test).*(1-interpFrac_onemat(:)),'omitnan')./normalizationfactor;
            end
        end
    end
end

[~,imax]=max(radonout(:));
[imaxpar,imaxperp,imaxphirot,imaxtshift]=ind2sub(size(radonout),imax);
twavefrontplot=squeeze(twavefront(:,:,imaxpar,imaxperp,imaxphirot))+tshiftvals(imaxtshift);
twavefrontplot(~latmap)=NaN; % only plot on considered lateral ranges

outRadonEllipse.cPar=cParvals(imaxpar);
outRadonEllipse.cPerp=cPerpvals(imaxperp);
outRadonEllipse.phirot=phirotvals(imaxphirot);
outRadonEllipse.tshift=tshiftvals(imaxtshift);
outRadonEllipse.cvals=squeeze(cvalsall(:,1,imaxpar,imaxperp,imaxphirot));
outRadonEllipse.twavefront=twavefrontplot;

if outRadonEllipse.cPar<outRadonEllipse.cPerp % if along>across, flip values
    if outRadonEllipse.phirot>=0
        outRadonEllipse.phirot=outRadonEllipse.phirot-90;
    else
        outRadonEllipse.phirot=outRadonEllipse.phirot+90;
    end
    tmphold=outRadonEllipse.cPar;
    outRadonEllipse.cPar=outRadonEllipse.cPerp;
    outRadonEllipse.cPerp=tmphold;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot3DRadonSum(out3DSWS,setupdataprocessing,data3D,data3Dfull,outplanes,setupSWS,setupMakePlanes,dataDir)

data3Dsurf(1:setupdataprocessing.nangles,:,:)=data3D;
data3Dsurf(setupdataprocessing.nangles+1,:,:)=data3D(1,:,:); % so full circle
data3Dsurffull(1:setupdataprocessing.nangles,:,:)=data3Dfull;
data3Dsurffull(setupdataprocessing.nangles+1,:,:)=data3Dfull(1,:,:); % so full circle
[rsurf,thetasurf]=meshgrid(outplanes(1).latmm,[setupdataprocessing.anglesDeg setupdataprocessing.anglesDeg(1)]);
xsurf=rsurf.*cosd(thetasurf);
ysurf=rsurf.*sind(thetasurf);
thetaradon=[setupdataprocessing.anglesDeg setupdataprocessing.anglesDeg(1)];

tmsstart=find(outplanes(1).tms>0,1,'first');

f2=figure(2);f2.Position= [715 43 1148 429];f2.Color=[1 1 1];
[savefolder,SWSsettingsname]=GenerateSaveFileName(dataDir,setupdataprocessing);
filenamegif=[savefolder 'Radon3DGif' SWSsettingsname '.gif'];
flagfirstframe=1;
for i=tmsstart:length(outplanes(1).tms)
    figure(2);subplot(1,2,1)
    img1=imagesc(setupdataprocessing.anglesDeg,outplanes(1).latmm,squeeze(data3Dsurffull(:,:,i))');hold on
    %img2= imagesc(setupdataprocessing.anglesDeg,outplanes(1).latmm,squeeze(data3Dsurf(:,:,i))','AlphaData',~isnan(squeeze(data3Dsurf(:,:,i))'));
    img1.AlphaData=0.5;

    tmp=gca;tmp.YDir='normal';
    hold on
    plot(setupdataprocessing.anglesDeg,out3DSWS.EachSWSEstimate.cvals.*outplanes(1).tms(i),'w-','LineWidth',2)
    title(['tms ' num2str(outplanes(1).tms(i))])
    caxis(setupMakePlanes.climsforplane)
    ylim(outplanes(1).latmm(end).*[0 1])
    hold off

    figure(2);subplot(1,2,2)
    s=surf(xsurf,ysurf,zeros(size(xsurf)),squeeze(data3Dsurf(:,:,i)));
    hold on
    s2=surf(xsurf,ysurf,zeros(size(xsurf)),squeeze(data3Dsurffull(:,:,i)));
    rradon=[out3DSWS.EachSWSEstimate.cvals'.*(outplanes(1).tms(i)-out3DSWS.EachSWSEstimate.tshift) out3DSWS.EachSWSEstimate.cvals(1).*(outplanes(1).tms(i)-out3DSWS.EachSWSEstimate.tshift)];
    xradon=rradon.*cosd(thetaradon);
    yradon=rradon.*sind(thetaradon);

    plot(xradon,yradon,'w','LineWidth',2)
    hold off

    caxis(setupMakePlanes.climsforplane)
    s.FaceColor='interp';
    s.EdgeColor = 'none';
    s2.FaceColor='interp';
    s2.EdgeColor = 'none';
    s2.FaceAlpha=.5;
    view(2)
    axis square
    xlim(outplanes(1).latmm(end).*[-1 1]);ylim(outplanes(1).latmm(end).*[-1 1])
    pause(.1);

    gifify(f2,flagfirstframe,filenamegif);flagfirstframe=2;
end

set(0,'CurrentFigure',1);clf;

nplot=2*setupdataprocessing.nangles;
ncol=round((floor(sqrt(nplot*2.5)))/2)*2;
nrow=ceil(nplot./ncol);

for iangle=1:length(outplanes)
    subplot(nrow,ncol,iangle*2-1);hold on
    imagesc(outplanes(iangle).tms,outplanes(iangle).latmm,outplanes(iangle).plane);
    caxis(setupMakePlanes.climsforplane)

    plot(out3DSWS.EachSWSEstimate.twavefront(iangle,:),outplanes(iangle).latmm,'w')
    title([num2str(out3DSWS.EachSWSEstimate.cvals(iangle),2) ' m/s'])
    tmp=gca;tmp.YDir='reverse';
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
polarplot(setupdataprocessing.anglesDeg*pi/180,out3DSWS.EachSWSEstimate.cvals,'k-') %plot SWS data from fit

polarplot([1 1]*out3DSWS.phiRot.*pi/180,[0 6],'k--')

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