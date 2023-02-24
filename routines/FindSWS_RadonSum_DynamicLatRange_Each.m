function [out,setupSWS]=FindSWS_RadonSum_DynamicLatRange_Each(outplanes,setupdataprocessing,setupSWS,setupMakePlane)
% Estimate SWS at each rotation angle using the dynamic lateral range Radon
% Sum method
% outplanes is a structure containing space-time planes for each acquistion rotation angle
% setupdataprocessing is a structure containing instructions for fitting
% setupSWS is a stucture containing information for fitting, and is modified and passed back out for record keeping
% setupMakePlane is a structure with plotting information
% out is a structure contains SWS estimates at each rotation angle

%% set Radon Params
setupSWS=setRadonSumParams(setupdataprocessing,setupSWS);

[setupSWS,ensemble_latrangestart,ensemble_latrangeend,latrangeheatmap_SH,latrangeheatmap_SecondSW,latrangeheatmap_ensemble]=SetupEnsemble(setupSWS,outplanes(1).latmm,setupdataprocessing.nangles);
%% Process each plane for SWS
for iangle=1:length(outplanes)
    [out(iangle)] = DynamicRangeRadonSum(outplanes(iangle),setupSWS,ensemble_latrangestart,ensemble_latrangeend);
end

%% Plot

if setupSWS.fignum

    [latrangeheatmap_SH,latrangeheatmap_SecondSW]=GenerateHeatMaps(out,latrangeheatmap_SH,latrangeheatmap_SecondSW);

    set(0,'CurrentFigure',1);clf;

    nplot=2*setupdataprocessing.nangles;
    ncol=round((floor(sqrt(nplot*2.5)))/2)*2;
    nrow=ceil(nplot./ncol);
    for iangle=1:length(out)
        figure(1);subplot(nrow,ncol,iangle*2-1)
        imagesc(outplanes(iangle).tms,outplanes(iangle).latmm,outplanes(iangle).plane);caxis(setupMakePlane.climsforplane)
        hold on;
        for i=1:length(out(iangle).ensemblewaves.tms_start)
            plot([out(iangle).ensemblewaves.tms_start(i),out(iangle).ensemblewaves.tms_end(i)],[outplanes(iangle).latmm(1) outplanes(iangle).latmm(end)],'r-')
        end
        tmp=gca;tmp.XTick=[];tmp.YTick=[];

        titlestring=[num2str(out(iangle).SH.speed,3) '\pm ' num2str(out(iangle).SH.speed_std,3)];
        title(titlestring);clear titlestring;

        %make plot of lat ranges used
        subplot(nrow,ncol,iangle*2)
        plot(latrangeheatmap_SH(iangle,:),outplanes(iangle).latmm)
        tmp=gca;tmp.YAxis.Direction='reverse';tmp.XTick=[];tmp.YTick=[];

    end

    figure(2)
    for i=1:setupdataprocessing.nangles;speedsSH(i)=out(i).SH.speed;errorbarsSH(i)=out(i).SH.speed_std;end
    for i=1:setupdataprocessing.nangles;speedsSecondSW(i)=out(i).SecondSW.speed(1);errorbarsSecondSW(i)=out(i).SecondSW.speed_std(1);end

    subplot(4,1,1)
    imagesc(1:setupdataprocessing.nangles,outplanes(1).latmm,(latrangeheatmap_SH./latrangeheatmap_ensemble')')
    colorbar;
    ylabel('latmm');

    subplot(4,1,2)
    errorbar(1:setupdataprocessing.nangles,speedsSH,errorbarsSH);
    ylabel('SH Speed')

    subplot(4,1,3)
    imagesc(1:setupdataprocessing.nangles,outplanes(1).latmm,(latrangeheatmap_SecondSW./latrangeheatmap_ensemble')')
    colorbar;
    ylabel('latmm')

    subplot(4,1,4)
    errorbar(1:setupdataprocessing.nangles,speedsSecondSW,errorbarsSecondSW);
    ylabel('SecondSWS speed')

    for i=1:4;subplot(4,1,i);tmp=gca;tmp.Position(3)=.775;end

end

end

%% ==========================================================================

function  setupSWS=setRadonSumParams(setupdataprocessing,setupSWS)

valstring=regexp(setupdataprocessing.SWSEstimationParams,'_originthres(\d+\.?\d?)minlatrange(\d+\.?\d?).*','tokens');
originthresval=str2double(valstring{1}{1});
minlatrangeval=str2double(valstring{1}{2});

% SWS params
setupSWS.originthresval=originthresval;

%setupSWS.minspeed; % already set
%setupSWS.maxspeed; % already set
setupSWS.tmsminstarttime=0;
setupSWS.tmsmaxstartandend=[100 100];

% additional params for ensemble
setupSWS.ensemble.latrangelength=2:10;
setupSWS.ensemble.maxlat=18;
setupSWS.ensemble.latrangestart=2:1:(setupSWS.ensemble.maxlat-min(setupSWS.ensemble.latrangelength));

setupSWS.ensemble.clustermaxlatstart=8;
setupSWS.ensemble.clusterminlatrange=minlatrangeval;

end

%=====================================================+%

function [setupSWS,ensemble_latrangestart,ensemble_latrangeend,latrangehotmap_SH,latrangehotmap_SecondSW,latrangehotmap_ensemble]=SetupEnsemble(setupSWS,latmm,nangles)

% create ensemble to check
ensemble_latrangestart=[];
ensemble_latrangeend=[];

for irangelength=1:length(setupSWS.ensemble.latrangelength)
    for istart=1:length(setupSWS.ensemble.latrangestart)
        if setupSWS.ensemble.latrangestart(istart)+setupSWS.ensemble.latrangelength(irangelength)<=setupSWS.ensemble.maxlat

            ensemble_latrangestart(length(ensemble_latrangestart)+1)=setupSWS.ensemble.latrangestart(istart);
            ensemble_latrangeend(length(ensemble_latrangeend)+1)=setupSWS.ensemble.latrangestart(istart)+setupSWS.ensemble.latrangelength(irangelength);
        end
    end
end
ensemble_latrangestart=round(ensemble_latrangestart,4);
ensemble_latrangeend=round(ensemble_latrangeend,4);
latrangehotmap_ensemble=zeros(length(latmm),1);
for i=1:length(ensemble_latrangestart)
    ilatstart = find(latmm>=ensemble_latrangestart(i),1,'first');
    ilatend   = find(latmm<=ensemble_latrangeend(i),1,'last');
    latrangehotmap_ensemble(ilatstart:ilatend)=latrangehotmap_ensemble(ilatstart:ilatend)+1;
end

latrangehotmap_SH=zeros(nangles,length(latmm));
latrangehotmap_SecondSW=zeros(nangles,length(latmm));

end

%===========================================================%

function [out] = DynamicRangeRadonSum(outplanes,setupSWS,ensemble_latrangestart,ensemble_latrangeend)

% run SWS fits
for ilatrange=1:length(ensemble_latrangestart)
    [~,outensemble(ilatrange).allWaves] = Calc_RadonSum_ForDynamicLatRange(outplanes.plane,outplanes.latmm,outplanes.tms,setupSWS,ensemble_latrangestart(ilatrange),ensemble_latrangeend(ilatrange),0,-5,5);
end
out.ensemblewaves=DetermineSpeedFromEnsemble(outensemble,outplanes,setupSWS.ensemble.clustermaxlatstart,setupSWS.ensemble.clusterminlatrange);
[out.SH,out.SecondSW]=AssignMultiwaves(out.ensemblewaves);

end

%===========================================================%

function out=DetermineSpeedFromEnsemble(outensemble,outplane,maxlatstart,minlatrange)

maxnumberofwaves=1;
tms_start_all=NaN(length(outensemble),maxnumberofwaves);
tms_end_all=NaN(length(outensemble),maxnumberofwaves);
speed_all=NaN(length(outensemble),maxnumberofwaves);
tms_start_all_shiftat0=NaN(length(outensemble),maxnumberofwaves);
tms_end_all_shiftatmaxlat=NaN(length(outensemble),maxnumberofwaves);
ilatstart=NaN(length(outensemble),maxnumberofwaves);
ilatend=NaN(length(outensemble),maxnumberofwaves);
datapeak=NaN(length(outensemble),maxnumberofwaves);


for i=1:length(outensemble)
    iwaves=1:min(maxnumberofwaves,length(outensemble(i).allWaves.datapeak));

    tms_start_all(i,iwaves)=outensemble(i).allWaves.tms_start(iwaves);
    tms_end_all(i,iwaves)=outensemble(i).allWaves.tms_end(iwaves);
    speed_all(i,iwaves)=outensemble(i).allWaves.speed(iwaves);


    tms_start_all_shiftat0(i,iwaves)=outensemble(i).allWaves.tms_start(iwaves)-outensemble(i).allWaves.latmm([outensemble(i).allWaves.ilatstart])./outensemble(i).allWaves.speed(iwaves);
    tms_end_all_shiftatmaxlat(i,iwaves)=outensemble(i).allWaves.tms_end(iwaves)+(outensemble(i).allWaves.latmm(end)-outensemble(i).allWaves.latmm([outensemble(i).allWaves.ilatend]))./outensemble(i).allWaves.speed(iwaves);

    ilatstart(i,iwaves)=[outensemble(i).allWaves.ilatstart];
    ilatend(i,iwaves)=[outensemble(i).allWaves.ilatend];

    datapeak(i,iwaves)=outensemble(i).allWaves.datapeak(iwaves);
end


% determine how many lateral positions a given waves is 'valid' for (has
% % similar wave'

speed_all_good=speed_all;
valididx=~isnan(speed_all_good);
tms_start_all=tms_start_all(valididx); %keep these for plotting
tms_end_all=tms_end_all(valididx);
ilatstart=ilatstart(valididx);ilatend=ilatend(valididx);

X=[tms_start_all_shiftat0(valididx),tms_end_all_shiftatmaxlat(valididx),speed_all_good(valididx)];

X(:,4)=(X(:,1)+X(:,2))/2;

numptsrequire=2; % threshold for number of close neighbors each point must have
numptsrequireincluster=floor(length(outensemble)*1/4);

clusterer = clusterDBSCAN('MinNumPoints',numptsrequire,'Epsilon',[.1 .5], ...
    'EnableDisambiguation',false);
[idx2,~] = clusterer(X(:,[3 4])); % speed and mean

[idx2]=ImposeLatRangeMinimumOnClusters(idx2,outplane.latmm(ilatstart),outplane.latmm(ilatend),maxlatstart,minlatrange,numptsrequireincluster);

%PlotClustering(X,idx2,outplane,tms_start_all,tms_end_all,ilatstart,ilatend); % allows to visualize clustering

for i=1:max(idx2) %only look at positive categories
    tms_start_ensemble(i)=mean(X(idx2==i,1));
    tms_start_ensemble_std(i)=std(X(idx2==i,1));

    tms_end_ensemble(i)=mean(X(idx2==i,2));
    tms_end_ensemble_std(i)=std(X(idx2==i,2));

    speed_ensemble(i)=mean(X(idx2==i,3));
    speed_ensemble_std(i)=std(X(idx2==i,3));

    numIDsincluded(i)=sum(idx2==i);

    ilatrange{i}(1,:)=ilatstart(idx2==i);
    ilatrange{i}(2,:)=ilatend(idx2==i);
end

%alternate way to consider speed
%speed_ensemble=(max(outensemble(1).allWaves.latmm))./(tms_end_ensemble-tms_start_ensemble)

%incase no groups detected
if max(idx2)<1
    tms_end_ensemble=NaN;tms_start_ensemble_std=NaN;
    tms_start_ensemble=NaN; tms_end_ensemble_std=NaN;
    speed_ensemble=NaN;speed_ensemble_std=NaN;
    numIDsincluded=NaN;
    ilatrange{1}=[NaN NaN]';
end


% always output from last ending wave to first ending wave -- approx SH first
[~,iorderwaves]=sort(tms_end_ensemble,'descend');

out.speed=speed_ensemble(iorderwaves);
out.speed_ensemble_std=speed_ensemble_std(iorderwaves);
out.tms_start=tms_start_ensemble(iorderwaves);
out.tms_start_std=tms_start_ensemble_std(iorderwaves);
out.tms_end=tms_end_ensemble(iorderwaves);
out.tms_end_std=tms_end_ensemble_std(iorderwaves);
out.numIDsincluded=numIDsincluded(iorderwaves);
for i=1:length(iorderwaves)
    out.ilatrange{i}=ilatrange{iorderwaves(i)};
end


end

function  [SH,SecondSW]=AssignMultiwaves(ensemblein)
SH.speed=ensemblein.speed(1);
SH.speed_std=ensemblein.speed_ensemble_std(1);

SH.tms_start=ensemblein.tms_start(1);
SH.tms_start_std=ensemblein.tms_start_std(1);

SH.tms_end=ensemblein.tms_end(1);
SH.tms_end_std=ensemblein.tms_end_std(1);

SH.numIDsincluded=ensemblein.numIDsincluded(1);
SH.ilatrange=ensemblein.ilatrange{1};
SH.hittingedgeflag=0;
SH.qualmetric=SH.speed_std;

if length(ensemblein.speed)==1
    SecondSW.speed=NaN;SecondSW.speed_std=NaN;
    SecondSW.tms_start=NaN;SecondSW.tms_start_std=NaN;
    SecondSW.tms_end=NaN;SecondSW.tms_end_std=NaN;
    SecondSW.numIDsincluded=NaN;SecondSW.ilatrange={};
else

    iSSW=2:length(ensemblein.speed);
    SecondSW.speed=ensemblein.speed(iSSW);
    SecondSW.speed_std=ensemblein.speed_ensemble_std(iSSW);

    SecondSW.tms_start=ensemblein.tms_start(iSSW);
    SecondSW.tms_start_std=ensemblein.tms_start_std(iSSW);

    SecondSW.tms_end=ensemblein.tms_end(iSSW);
    SecondSW.tms_end_std=ensemblein.tms_end_std(iSSW);

    SecondSW.numIDsincluded=ensemblein.numIDsincluded(iSSW);
    for i=1:length(iSSW)
        SecondSW.ilatrange{i}=ensemblein.ilatrange{i};
    end

end

end

%===========================================================%

function [newidx2]=ImposeLatRangeMinimumOnClusters(idx2,latstart,latend,maxlatstart,minlatrange,numminpoints_incluster)

categorieslabeloptions=unique(idx2);
for i1=1:length(categorieslabeloptions)
    iwaveincluster=find(idx2==categorieslabeloptions(i1));
    if min(latstart(iwaveincluster))>maxlatstart %if doesn't start early enough
        idx2(iwaveincluster)=-1; %remove cluster as viable option
    end
    if max(latend(iwaveincluster))-min(latstart(iwaveincluster))<minlatrange % if total range is too small
        idx2(iwaveincluster)=-1; %remove cluster as viable option
    end
    if length(iwaveincluster)<numminpoints_incluster % if overall cluster doesn't have enough points
        idx2(iwaveincluster)=-1; %remove cluster as viable option
    end
end

%move category labels down to lowest numbers
categorieslabeloptions=unique(idx2);
newidx2=idx2;
categorieslabeloptionsabovezero=categorieslabeloptions(categorieslabeloptions>0);

for i=1:length(categorieslabeloptionsabovezero)
    newidx2(idx2==categorieslabeloptionsabovezero(i))=i; % replace with new number, stepping up
end

end

%===========================================================%

function PlotClustering(X,idx2,outplane,tms_start_all,tms_end_all,ilatstart,ilatend)

categorieslabeloptions=unique(idx2);

if any(idx2==-1)
    colorlist=[0 0 0 ;0 1 0;1 0 0;0 0 1;.5 1 .5; 1 .5 .5];
else
    colorlist=[0 1 0;1 0 0;0 0 1;.5 1 .5; 1 .5 .5];
end
if size(colorlist,1)<length(categorieslabeloptions)
    colorlist=[colorlist; jet(length(categorieslabeloptions)-size(colorlist,1))];
end
colorlist=colorlist(1:length(categorieslabeloptions),:);

figure(2);clf;
subplot(3,2,1)
gscatter(X(:,1),X(:,2),idx2,colorlist);xlim([-10 45]);ylim([-10 45]);axis square
xlabel('start');ylabel('end')
subplot(3,2,2)
gscatter(X(:,3),X(:,4),idx2,colorlist);xlim([0 6]);ylim([0 20]);axis square
xlabel('speed');ylabel('mean')

subplot(3,2,3)
imagesc(outplane.tms,outplane.latmm,outplane.plane);axis image;hold on
clim([-.5 .5])
subplot(3,2,4)
imagesc(outplane.tms,outplane.latmm,outplane.plane);axis image;hold on
clim([-.5 .5])
subplot(3,2,6)
imagesc(outplane.tms,outplane.latmm,outplane.plane);axis image;hold on
clim([-.5 .5])

subplot(3,2,5);hold on;
colormaplat=jet(size(X,1));
for i=1:size(X,1);plot(X(i,3),X(i,4),'.','Color',colormaplat(i,:));end
xlim([-10 45]);ylim([-10 45])

figure(2);
for i1=1:length(categorieslabeloptions)
    iwaveincluster=find(idx2==categorieslabeloptions(i1));
    if categorieslabeloptions(i1)==-1
        subplot(3,2,4)
        for i2=1:length(iwaveincluster)
            plot([tms_start_all(iwaveincluster(i2)) tms_end_all(iwaveincluster(i2))],outplane.latmm([ilatstart(iwaveincluster(i2)) ilatend(iwaveincluster(i2))]),'-','Color',[.5 .5 .5])
        end
        subplot(3,2,6)
        for i2=1:length(iwaveincluster)
            plot([X(iwaveincluster(i2),1) X(iwaveincluster(i2),2)],[0 outplane.latmm(end)],'-','Color',colorlist(i1,:))
        end


    else
        subplot(3,2,3)
        for i2=1:length(iwaveincluster)
            plot([tms_start_all(iwaveincluster(i2)) tms_end_all(iwaveincluster(i2))],outplane.latmm([ilatstart(iwaveincluster(i2)) ilatend(iwaveincluster(i2))]),'-','Color',colorlist(i1,:))
        end

        subplot(3,2,6)
        for i2=1:length(iwaveincluster)
            plot([X(iwaveincluster(i2),1) X(iwaveincluster(i2),2)],[0 outplane.latmm(end)],'-','Color',colorlist(i1,:))
        end

        subplot(3,2,4)
        plot([mean(X(iwaveincluster,1)) mean(X(iwaveincluster,2))],[0 outplane.latmm(end)],'--','Color',colorlist(i1,:))
    end
end

end

function [latrangeheatmap_SH,latrangeheatmap_SecondSW]=GenerateHeatMaps(out,latrangeheatmap_SH,latrangeheatmap_SecondSW)

for iangle=1:length(out)
    for i=1:size(out(iangle).SH.ilatrange,2)
        if ~isnan(out(iangle).SH.ilatrange(1,i))
            tmp=out(iangle).SH.ilatrange(1,i):out(iangle).SH.ilatrange(2,i);
            latrangeheatmap_SH(iangle,tmp)=latrangeheatmap_SH(iangle,tmp)+1;%./length(tmp)
        end
    end
    if ~isempty(out(iangle).SecondSW.ilatrange)
        for i=1:size(out(iangle).SecondSW.ilatrange{1},2)
            tmp=out(iangle).SecondSW.ilatrange{1}(1,i):out(iangle).SecondSW.ilatrange{1}(2,i);
            latrangeheatmap_SecondSW(iangle,tmp)=latrangeheatmap_SecondSW(iangle,tmp)+1;%./length(tmp);
        end
    end

end
end