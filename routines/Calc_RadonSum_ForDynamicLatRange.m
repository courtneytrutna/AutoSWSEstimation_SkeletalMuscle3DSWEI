function [out,out2] = Calc_RadonSum_ForDynamicLatRange(plane,latmm,tms,setupSWS,minlat,maxlat,fignum,color_min, color_max)

% Calcuates Radon Sum over a single plane and estimates SWS
% out contains information from primary peak (global max), out2 contains information
% from additional local peaks in the Radon Sum space

% this script is very similar to Calc_RadonSum, with changes to input
% parameters for use in the dynamic range method

prominencethreshold=1;

tmsminstarttime=setupSWS.tmsminstarttime;
tmsmaxstartandend=setupSWS.tmsmaxstartandend;
origintmsthreshold=setupSWS.originthresval;
minspeed=setupSWS.minspeed;
maxspeed=setupSWS.maxspeed;

% Define range over which to latSum
ilatstart = find(latmm>=minlat,1,'first');
ilatend   = find(latmm<=maxlat,1,'last');

% Find time start- cut out pre-zero start points
[~, minstarttimeindex]=min(abs(tms-tmsminstarttime(1)));
[~, maxstarttimeindex]=min(abs(tms-tmsmaxstartandend(1)));
[~, maxendtimeindex]=min(abs(tms-tmsmaxstartandend(2)));

[latsums]=FindLatSumsUsingEndpoints(plane,ilatstart,ilatend,minstarttimeindex,maxstarttimeindex,maxendtimeindex);
latsumsoriginal=latsums;

[latsums,itstart,itend,tms_start,tms_end,hittingedgeflag,ilocaltstart,ilocaltend,tms_start_local,tms_end_local,latsums_sub] = FindBestAndLocalPeaks(latsums,tms,latmm,ilatstart,ilatend,minstarttimeindex,maxstarttimeindex,maxendtimeindex,prominencethreshold,minspeed,maxspeed,origintmsthreshold);

alphalatsums=ones(size(latsums))*.5;alphalatsums(~isnan(latsums))=0;

slope = (tms_end-tms_start)/(ilatend-ilatstart);
speed = mean(diff(latmm))/slope;


if fignum
    figure(fignum)
    clf
    subplot(1,2,1)
    imagesc(tms,latmm,plane)
    hold on
    plot([tms_start tms_end],[latmm(ilatstart) latmm(ilatend)],'w', 'LineWidth',2)
    for i=1:length(ilocaltstart);plot([tms_start_local(i) tms_end_local(i)],[latmm(ilatstart) latmm(ilatend)],'r');end

    hold off
    xlabel('time (ms)');
    ylabel('lateral position (mm)');
    title(['gSWS = ' num2str(speed) ' m/s'])
    caxis([color_min color_max])
    colorbar

    subplot(1,2,2)
    imagesc(tms(1:(end-1)),tms(1:(end-1)),latsums)
    hold on
    imagesc(tms(1:(end-1)),tms(1:(end-1)),latsumsoriginal,'AlphaData',alphalatsums)
    plot([tms_end_local tms_end_local],[tms_start_local tms_start_local],'.r','MarkerSize',12)
    plot([tms_end tms_end],[tms_start tms_start],'.w','MarkerSize',16)
    hold off
    set(gca,'DataAspectRatio',[1 1 1])
    xlabel('end time (ms)')
    ylabel('start time (ms)')
end

if ~isempty(ilocaltstart)
    datapeak = latsums_sub(1);
else
    datapeak = NaN;
end

out.datapeak  = datapeak;
out.speed     = speed;
out.itstart   = itstart;
out.itend     = itend;
out.ilatstart = ilatstart;
out.ilatend   = ilatend;
out.latsums   = latsums;
out.plane     = plane;
out.latmm     = latmm;
out.tms       = tms;
out.tms_start = tms_start;
out.tms_end = tms_end;
out.hittingedgeflag= hittingedgeflag;
out.ttpqualmetric= [];%ttpqualmetric;

out2.datapeak = latsums_sub;
out2.ilocaltstart = ilocaltstart;
out2.ilocaltend   = ilocaltend;
out2.ilatstart    = ilatstart;
out2.ilatend      = ilatend;
out2.latsums      = latsums;
out2.latmm        = latmm;
out2.tms          = tms;
out2.tms_start    = tms_start_local;
out2.tms_end      = tms_end_local;

slope = (out2.tms_end-out2.tms_start)/(out2.ilatend-out2.ilatstart);
out2.speed = mean(diff(latmm))./slope;
end

%==========================================================================

function [latsums]=FindLatSumsUsingEndpoints(plane,ilatstart,ilatend,minstarttimeindex,maxstarttimeindex,maxendtimeindex)

[nlats,ntimes]=size(plane);

interpIndx=zeros(ntimes-1,nlats);         % get indices and fractions for
interpFrac=zeros(ntimes-1,nlats);         %   interpolation at lat positions

ilats=ilatstart:ilatend;

for idiff=1:ntimes-1   %run across every time step

    tvals=(idiff-1)/(length(ilats)-1)*(ilats-ilats(1));  % 'ideal' time steps across each lat to connect the way we'd want to
    interpIndx(idiff,ilats)=floor(tvals);
    interpFrac(idiff,ilats)=1-(tvals-interpIndx(idiff,ilats));
end

latsums = nan(ntimes-1,ntimes-1);

maxistart=min(maxstarttimeindex,ntimes-1); %use given index, or last index available if smaller
maxiend=min(maxendtimeindex,ntimes-1);

startTimeIndex=minstarttimeindex;
endoffset = 0;

for istart=startTimeIndex:maxistart
    iendvector=istart+endoffset:maxiend;
    idiff=iendvector-istart+1;            % +1 for matlab numbering  Find difference between beginning and end time steps

    idx0vals=interpIndx(idiff,ilats)+istart;  %Index zero vals
    fracs=interpFrac(idiff,ilats);

    indices1 = ones(size(iendvector))'*ilats + (idx0vals-1)*size(plane,1); %Each 1x28 matrices of indexes-- drawing the 'line'
    indices2= ones(size(iendvector))'*ilats+ (idx0vals)*size(plane,1); % (id0vals+1)-1

    % Above is faster than using sub2ind
    % indices2test = sub2ind(size(plane),ilats,idx0vals+1);
    %indices1test= sub2ind(size(plane),ilats,idx0vals);

    latsums(istart,iendvector) = sum(plane(indices1).*fracs+plane(indices2).*(1-fracs),2,'omitnan')';  % sum over lat pos,weighting the data points on either side of the 'line' appropriately
end

end

%==========================================================================

function latsums = LimitMaxSpeed(latsums,t,lat,ilatstart,ilatend,minspeed,maxspeed)

% input lat in mm, t in ms

%    [nlats,ntimes]=size(plane);
%
%    ilatstart=1;
%    ilatend=nlats;

dlat = mean(diff(lat));
dt   = mean(diff(t));

% fastest could travel, in "indicies" (ie if max speed is 1 full distance/3 time steps, remove lines that were only two time steps
deltait_max = floor(dlat*(ilatend-ilatstart)/(maxspeed*dt));    % dt in ms
deltait_min = floor(dlat*(ilatend-ilatstart)/(minspeed*dt));    % dt in ms

[nsums1,nsums2]=size(latsums);
if nsums1~=nsums2
    error('error')
end

col=(1:nsums1)';
itstartmat=repmat(col,1,nsums2);
itendmat=itstartmat';

latsums(itendmat<=itstartmat+deltait_max)=NaN;  % replace too high speeds with NaN, including negative speeds
latsums(itendmat>=itstartmat+deltait_min)=NaN;  % replace too high speeds with NaN, including negative speeds
end

%==========================================================================

function latsums = LimitBackProjectionToOrigin(latsums,tms,latmm,ilatstart,ilatend,origintmsthreshold)
% input lat in mm, t in ms

lineequation_tmsatlat0=-(tms'-tms)./(latmm(ilatend)-latmm(ilatstart))*latmm(ilatstart)+tms;

mask=abs(lineequation_tmsatlat0(2:end,2:end))>origintmsthreshold;

latsums(mask')=NaN; %Nan out lines that dont fit

end

%==========================================================================

function [latsums,itstart,itend,tms_start,tms_end,hittingedgeflag,ilocaltstart_sub,ilocaltend_sub,tms_start_local,tms_end_local,latsums_sub]=FindBestAndLocalPeaks(latsums,tms,latmm,ilatstart,ilatend,minstarttimeindex,maxstarttimeindex,maxendtimeindex,prominencethreshold,minspeed,maxspeed,origintmsthreshold)

% multipeaks
latsumstmp=latsums;
latsumstmp(isnan(latsums))=min(latsums(:));
[i1,p1]=islocalmax(latsumstmp,1);[i2,p2]=islocalmax(latsumstmp,2);
flaglocalpeaks=and(i1,i2);
ilocalpeaks=find(flaglocalpeaks);
vallocalpeaks=latsums(ilocalpeaks);
[~,iorder]=sort(vallocalpeaks,'descend');
ilocalpeaks=ilocalpeaks(iorder);
[ilocaltstart, ilocaltend]=ind2sub(size(latsums),ilocalpeaks);

% was creating problem with prominence and late lateral positions-- small
% ROI because of projection back to origin skewing the promienence calc
% solution: find multipeaks and prominence, *then* limit with speed,
% project to origin. Below will remove because in NaN territory

latsums = LimitMaxSpeed(latsums,tms,latmm,ilatstart,ilatend,minspeed,maxspeed);
latsums = LimitBackProjectionToOrigin(latsums,tms,latmm,ilatstart,ilatend,origintmsthreshold);

% Remove 'bad' peask from multipeaks
idxremove=[];
for itest=1:length(ilocaltstart)
    % remove if hitting edge
    if ilocaltend(itest)>1 % need so -1 index checks work below
        % remove anything at edge of ROI
        if any([ilocaltstart(itest)==minstarttimeindex,ilocaltstart(itest)==maxstarttimeindex,ilocaltend(itest)==maxendtimeindex,...
                isnan(latsums(ilocaltstart(itest)+1,ilocaltend(itest))),isnan(latsums(ilocaltstart(itest)-1,ilocaltend(itest))),...
                isnan(latsums(ilocaltstart(itest),ilocaltend(itest)+1)),isnan(latsums(ilocaltstart(itest),ilocaltend(itest)-1)),...
                isnan(latsums(ilocaltstart(itest)+1,ilocaltend(itest)+1)),isnan(latsums(ilocaltstart(itest)-1,ilocaltend(itest)-1)),...
                isnan(latsums(ilocaltstart(itest)+1,ilocaltend(itest)-1)),isnan(latsums(ilocaltstart(itest)-1,ilocaltend(itest)+1))])
            idxremove(length(idxremove)+1)=itest;
        end
    else
        idxremove(length(idxremove)+1)=itest;
    end

    %remove sum-squared prominence less than 0.25
    if sqrt(p2(ilocaltstart(itest),ilocaltend(itest)).^2+p1(ilocaltstart(itest),ilocaltend(itest)).^2)<prominencethreshold
        idxremove(length(idxremove)+1)=itest;
    end
end
ilocaltstart(idxremove)=[];
ilocaltend(idxremove)=[];

idxremove=[];
for itest=1:length(ilocaltstart)
    %remove peaks that are too close to peak of high magnitude (particularly useful for digonals)
    % we probaly do want to remove peaks that are close to ROI-edge
    % peaks, so leave in this loop
    if any((sqrt((ilocaltstart(1:(itest-1))-ilocaltstart(itest)).^2)+sqrt((ilocaltend(1:(itest-1))-ilocaltend(itest)).^2))<3) % when itest=1, array will be empty, any([])=0
        idxremove(length(idxremove)+1)=itest;
    end
end
ilocaltstart(idxremove)=[];
ilocaltend(idxremove)=[];

idxremove=[];
for itest=1:length(ilocaltstart) % need to do after remove edge of ROI peaks above, because they cross alot

    % remove trajectories that 'cross' previous/higher trajectories
    % over lateral range tested
    % if cross, start-start and end-end will be different signs--
    % multiply and check for <0
    if any(((ilocaltstart(1:(itest-1))-ilocaltstart(itest)).*(ilocaltend(1:(itest-1))-ilocaltend(itest)))<0)
        idxremove(length(idxremove)+1)=itest;
    end

end
ilocaltstart(idxremove)=[];
ilocaltend(idxremove)=[];

%subsample around peaks
ww=2; % window for subsampling
latsums_sub=[];
ilocaltstart_sub=[];
ilocaltend_sub=[];
for isub=1:length(ilocaltstart)
    if and((ilocaltstart(isub)-ww)>=1,(ilocaltstart(isub)+ww)<=size(latsums,2))
        p=polyfit((ilocaltstart(isub)-ww):(ilocaltstart(isub)+ww),latsums((ilocaltstart(isub)-ww):(ilocaltstart(isub)+ww),ilocaltend(isub)),2);
        if ~isnan(p(1))
            ilocaltstart_sub(isub) = -p(2)/(2*p(1));
            ilocaltstart_sub_latsums(isub)=p(1)*(ilocaltstart_sub(isub))^2+p(2)*ilocaltstart_sub(isub)+p(3);
        else %hitting edge of ROI
            ilocaltstart_sub(isub)=ilocaltstart(isub);
            ilocaltstart_sub_latsums(isub)=latsums(ilocaltstart(isub),ilocaltend(isub));
        end
    else
        ilocaltstart_sub(isub)=ilocaltstart(isub);
        ilocaltstart_sub_latsums(isub)=latsums(ilocaltstart(isub),ilocaltend(isub));
    end

    if and((ilocaltend(isub)-ww)>=1,(ilocaltend(isub)+ww)<=size(latsums,2))
        p=polyfit((ilocaltend(isub)-ww):(ilocaltend(isub)+ww),latsums(ilocaltstart(isub),(ilocaltend(isub)-ww):(ilocaltend(isub)+ww)),2);
        ilocaltend_sub(isub) = -p(2)/(2*p(1));
        ilocaltend_sub_latsum(isub)=p(1)*(ilocaltend_sub(isub))^2+p(2)*ilocaltend_sub(isub)+p(3);
    else
        ilocaltend_sub(isub)=ilocaltend(isub);
        ilocaltend_sub_latsum(isub)=latsums(ilocaltstart(isub),ilocaltend(isub));
    end

    latsums_sub(isub)=(ilocaltstart_sub_latsums(isub)+ilocaltend_sub_latsum(isub))/2; %just average to avoid having to do 2D fitting

end

if ~isempty(ilocaltstart)
    itstart=ilocaltstart_sub(1); % pass out highest peak as primary
    itend=ilocaltend_sub(1);
    hittingedgeflag=0;
else
    itstart=1;
    itend=1;
    hittingedgeflag=1;
end

tms_start=interp1(1:length(tms),tms,itstart);
tms_end=interp1(1:length(tms),tms,itend);

tms_start_local    = interp1(1:length(tms),tms,ilocaltstart);
tms_end_local      = interp1(1:length(tms),tms,ilocaltend);
end



