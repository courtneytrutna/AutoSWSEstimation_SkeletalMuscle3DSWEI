function [TTP,setupSWS]=Calc_TTP(outplanes,setupSWS)
% Calculate time-to-peak (TTP) for each lateral position for each plane
% outplanes and setupSWS are both structures with certain fields expected,
% see below
%
% TTP contains all time-to-peak datapoints,
% size of TTP: number of acquistion angles x number of lateral positions
%
% setupSWS is updated with TTP parameters (set below) and passed back out
% for record keeping


%% Set options
setupSWS.TTPparams.TTPrange=[.5 25];
setupSWS.TTPparams.TTPmaxspeed=setupSWS.maxspeed*1.1;
setupSWS.TTPparams.TTPminspeed=setupSWS.minspeed*.9;
setupSWS.TTPparams.TTPmedianfiltblocksize=5;

%calculate dt for later subsample step
dt=median(diff(outplanes(1).tms));

% mask region over which to look for peak
ttpmask=ones(size(outplanes(1).plane));
ttpmask(abs(outplanes(1).latmm'./outplanes(1).tms)>setupSWS.maxspeed*1.1)=0;
ttpmask(abs(outplanes(1).latmm'./outplanes(1).tms)<setupSWS.minspeed*0.9)=0;
ttpmask(ones(size(outplanes(1).latmm,2),1).*outplanes(1).tms<setupSWS.TTPparams.TTPrange(1))=0;
ttpmask(ones(size(outplanes(1).latmm,2),1).*outplanes(1).tms>setupSWS.TTPparams.TTPrange(2))=0;


%loop over all collected angles
for iangle=1:length(outplanes)

    [~,imax]=max(outplanes(iangle).plane.*ttpmask,[],2);

    % subsample around peak
    for ilat=1:length(imax)
        if or((imax(ilat)+1)>size(outplanes(iangle).plane,2),(imax(ilat)-1)<1)
            peakloc=0; %if at edge of plane, don't interpolate
        else %otherwise, do quadratic interpolation
            alpha=outplanes(iangle).plane(ilat,imax(ilat)-1);
            beta=outplanes(iangle).plane(ilat,imax(ilat));
            gamma=outplanes(iangle).plane(ilat,imax(ilat)+1);

            peakloc=1/2*(alpha-gamma)./(alpha-2*beta+gamma);
        end

        TTP(iangle,ilat)=outplanes(iangle).tms(imax(ilat))+peakloc*dt;
    end

end

%% Clean up TTP data

[latmesh]=repmat(outplanes(1).latmm,[size(TTP,1),1]);
TTP((TTP*setupSWS.TTPparams.TTPmaxspeed)<latmesh)=NaN; % remove unrealistically fast speeds
TTP((TTP*setupSWS.TTPparams.TTPminspeed)>latmesh)=NaN; % remove unrealistically slow speeds
TTP=medfilt1(TTP,setupSWS.TTPparams.TTPmedianfiltblocksize,[],2,'includenan'); % median filter

end