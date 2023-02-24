function [XCorrout,setupSWS]=Calc_XCorr(outplanes,setupSWS,setupdataprocessing)
% Calculate cross correlation time shifts for each lateral position for each plane
% outplanes, setupSWS, and setupdataprocessing are structures with certain fields expected,
% see below
%
% XCorrout contains all estimated timeshift datapoints,
% size of XCorrout: number of acquistion angles x number of lateral positions
%
% setupSWS is updated with XCorr parameters (set below) and passed back out
% for record keeping

%% Set options

% setup parameters
setupSWS=setXCorrParams(setupdataprocessing,setupSWS);
setupSWS.XCorrparams.XCorrtimerange=[.5 25];

%% calculate useful params
dt=median(diff(outplanes(1).tms));
ilatstep=round(setupSWS.XCorrparams.latstepsize_mm./median(diff(outplanes(1).latmm))); % needs to be an integer

if mod(ilatstep,2) %want to be even so can split in two -- if odd
    ilatstep=ilatstep+1; %make even
end

ilatstep_eachway=ilatstep/2; % this should be an integer and if it's not it should error
% update to true value
setupSWS.XCorrparams.latstepsize_mm=ilatstep*median(diff(outplanes(1).latmm));


%find time range to use for analysis
% cut out the first couple time steps to avoid noise influencing
% measurements
xcorrmask=ones(size(outplanes(1).plane));
xcorrmask(abs(outplanes(1).latmm'./outplanes(1).tms)>setupSWS.maxspeed*1.1)=0;
xcorrmask(abs(outplanes(1).latmm'./outplanes(1).tms)<setupSWS.minspeed*0.9)=0;
xcorrmask(ones(size(outplanes(1).latmm,2),1).*outplanes(1).tms<setupSWS.XCorrparams.XCorrtimerange(1))=0;
xcorrmask(ones(size(outplanes(1).latmm,2),1).*outplanes(1).tms>setupSWS.XCorrparams.XCorrtimerange(2))=0;

minshift=1./(setupSWS.maxspeed*1.1);
maxshift=1./(setupSWS.minspeed*.9);

%% loop over all collected angles
for iangle=1:length(outplanes)
    planetmp=outplanes(iangle).plane.*xcorrmask; % mask out things we don't want

    % run cross corr for all lat combos
    ilattest=(ilatstep_eachway+1):(length(outplanes(iangle).latmm)-ilatstep_eachway);
    % ilat=1 is all 0, so skip that one

    for ilat=ilattest
        % estimate shift
        %chose two lat positions for comparison
        tmp1=planetmp(ilat-ilatstep_eachway,:);
        tmp2=planetmp(ilat+ilatstep_eachway,:);
        % institute minimum signal threshold
        tmp1(tmp1<setupSWS.XCorrparams.signalthreshold)=setupSWS.XCorrparams.signalthreshold; tmp2(tmp2<setupSWS.XCorrparams.signalthreshold)=setupSWS.XCorrparams.signalthreshold;
        % calc cross corr
        [c,lags]=xcorr_fast(tmp2',tmp1');
        c(lags<0)=0; % remove negative moving waves
        [maxc,imax]=max(c);

        if imax==1
            shift_ms(ilat)=NaN;
            maxc_all(iangle,ilat)=0;
            continue
        end

        % subsample around peak
        try
            alpha=c(imax-1);
            beta=c(imax);
            gamma=c(imax+1);
            peakloc=1/2*(alpha-gamma)./(alpha-2*beta+gamma);
        catch
            peakloc=0;
        end

        % record shift, divided over mm distance of shift
        shift_ms(ilat)=(lags(imax)*dt+peakloc*dt)./setupSWS.XCorrparams.latstepsize_mm; % shift in ms, normalized by length comparing (in mm)
        maxc_all(iangle,ilat)=maxc;

    end
    XCorr(iangle,1:ilattest(1))=NaN;  % to keep sizes the same
    XCorr(iangle,ilattest)=shift_ms(ilattest); % note shift_ms is at ilattest, even if ilattest doesn't start at 0
    XCorr(iangle,ilattest(end):length(outplanes(1).latmm))=NaN; % to keep sizes the same
    clear shift_ms

end

maxc_all(iangle,1:ilattest(1))=NaN;  % to keep sizes the same
maxc_all(iangle,ilattest(end):length(outplanes(1).latmm))=NaN; % to keep sizes the same

maxc_all(XCorr<minshift)=NaN;
maxc_all(XCorr>maxshift)=NaN;
XCorr(XCorr<minshift)=NaN;
XCorr(XCorr>maxshift)=NaN;

XCorrout.XCorrshifts_mspermm=XCorr;%already normed
XCorrout.CorrCoeff=maxc_all;
end


function [c,lags] = xcorr_fast(x, y)
maxlag = max(size(x,1), size(y,1)) - 1;
c1 = crosscorr(x, y, maxlag);
c = c1 ./ sqrt(sum(abs(x).^2)*sum(abs(y).^2));
lags = -maxlag:maxlag;
end

function c = crosscorr(x,y,maxlag)
m = max(size(x,1), size(y,1));
mxl = min(maxlag, m-1);
m2 = 2^nextpow2(2*m-1);

X = fft(x, m2, 1);
Y = fft(y, m2, 1);
c1 = real(ifft(X.*conj(Y),[],1));
c = [c1(m2 - mxl + (1:mxl)); c1(1:mxl+1)];
end

function [setupSWS]=setXCorrParams(setupdataprocessing,setupSWS)

valstring=regexp(setupdataprocessing.SWSEstimationParams,'.*latstepmm(\d+\.?\d?)sigthreshold(-*\d+\.?\d?).*','tokens');
latstepsize_mm=str2double(valstring{1}{1});
signalthreshold=str2double(valstring{1}{2});
setupSWS.XCorrparams.latstepsize_mm=latstepsize_mm;
setupSWS.XCorrparams.signalthreshold=signalthreshold;
end
