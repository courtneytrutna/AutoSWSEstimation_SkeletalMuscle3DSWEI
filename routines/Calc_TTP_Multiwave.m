function [TTP,setupSWS]=Calc_TTP_Multiwave(outplanes,setupSWS)
% Calculate multiwave time-to-peak (TTP) for each lateral position for each plane
% outplanes and setupSWS are both structures with certain fields expected,
% see below
%
% this script is similar to Calc_TTP except times to the two largest local peaks
% are output, rather than just the global peak, allowing for multiwave
% detection
%
% TTP contains all time-to-peak datapoints,
% size of TTP: number of acquistion angles x number of lateral positions x number of peaks (2)
%
% setupSWS is updated with TTP parameters (set below) and passed back out
% for record keeping

warning('off','signal:findpeaks:largeMinPeakHeight')
%% Set options
setupSWS.TTPparams.TTPrange=[.5 25];
setupSWS.TTPparams.TTPmaxspeed=setupSWS.maxspeed*1.1;
setupSWS.TTPparams.TTPminspeed=setupSWS.minspeed*.9;
setupSWS.TTPparams.TTPmedianfiltblocksize=5;

%calculate dt for later subsample step
dt=median(diff(outplanes(1).tms));

% mask region over which to look for peak

ttpmask=ones(size(outplanes(1).plane));
ttpmask(abs(outplanes(1).latmm'./outplanes(1).tms)>setupSWS.TTPparams.TTPmaxspeed)=0;
ttpmask(abs(outplanes(1).latmm'./outplanes(1).tms)<setupSWS.TTPparams.TTPminspeed)=0;
ttpmask(ones(size(outplanes(1).latmm,2),1).*outplanes(1).tms<setupSWS.TTPparams.TTPrange(1))=0;
ttpmask(ones(size(outplanes(1).latmm,2),1).*outplanes(1).tms>setupSWS.TTPparams.TTPrange(2))=0;

TTP=NaN(length(outplanes),size(outplanes(1).plane,1),2); % up to two peaks
%loop over all collected angles
for iangle=1:length(outplanes)

    maskedplane=outplanes(iangle).plane.*ttpmask;

    % subsample around peak
    for ilat=1:size(maskedplane,1)
        [peaks,tmp1,~,p]=findpeaks(maskedplane(ilat,:),'NPeaks',2,'MinPeakProminence',.05,'MinPeakHeight',0.25,'SortStr','descend');

        if length(p)==2
            if p(1)*.25>p(2) % require prominence of second point to be at least a quarter the first
                tmp1(2)=[];
                peaks(2)=[];
            end
        end
        imax(ilat,1:length(tmp1))=tmp1;

        for i=1:length(peaks)
            if or((imax(ilat,i)+1)>size(outplanes(iangle).plane,2),(imax(ilat,i)-1)<1)
                peakloc=0; %if at edge of plane, don't interpolate
            else %otherwise, do quadratic interpolation
                alpha=outplanes(iangle).plane(ilat,imax(ilat,i)-1);
                beta=outplanes(iangle).plane(ilat,imax(ilat,i));
                gamma=outplanes(iangle).plane(ilat,imax(ilat,i)+1);

                peakloc=1/2*(alpha-gamma)./(alpha-2*beta+gamma);
            end

            TTP(iangle,ilat,i)=outplanes(iangle).tms(imax(ilat,i))+peakloc*dt;
        end
    end

end

%% Clean up TTP data

[latmesh]=repmat(outplanes(1).latmm,[size(TTP,1),1,size(TTP,3)]);
TTP((TTP*setupSWS.TTPparams.TTPmaxspeed)<latmesh)=NaN; % remove unrealistically fast speeds
TTP((TTP*setupSWS.TTPparams.TTPminspeed)>latmesh)=NaN; % remove unrealistically slow speeds
TTP=medfilt1(TTP,setupSWS.TTPparams.TTPmedianfiltblocksize,[],2,'includenan'); % median filter

end