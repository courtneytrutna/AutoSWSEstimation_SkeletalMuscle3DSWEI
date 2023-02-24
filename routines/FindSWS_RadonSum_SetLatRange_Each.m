function [out,setupSWS]=FindSWS_RadonSum_SetLatRange_Each(outplanes,setupdataprocessing,setupSWS,setupMakePlane,dataDir)

% Estimate SWS at each rotation angle using the Radon Sum method (fixed or angle informed lateral range)
% outplanes is a structure containing space-time planes for each acquistion rotation angle
% setupdataprocessing is a structure containing instructions for fitting
% setupSWS is a stucture containing information for fitting, and is modified and passed back out for record keeping
% setupMakePlane is a structure with plotting information
% dataDir is used to load B-mode information
% out is a structure contains SWS estimates at each rotation angle

%% set Radon Params
setupSWS=setRadonSumParams(setupdataprocessing,setupSWS,dataDir);

%% Process each plane for SWS
for iangle=1:length(outplanes)
    [out(iangle).SH,out(iangle).SecondSW] = Calc_RadonSum(outplanes(iangle).plane,outplanes(iangle).latmm,outplanes(iangle).tms,setupSWS.SH(iangle),0,setupMakePlane.climsforplane); %zero is fignum

    if regexp(setupdataprocessing.SWSEstimationParams,'multiwave')
        [out(iangle)]=CleanAndCategorizeMultiwaves(out(iangle),1);
    else
        [out(iangle)]=CleanAndCategorizeMultiwaves(out(iangle),0); %second zero prevents any SV SH flippping, but still clean up any multiple peaks found
    end
end
if setupSWS.fignum
    %% Plot
    figure(setupSWS.fignum)
    set(0,'CurrentFigure',1);clf;

    if setupdataprocessing.nangles<144
        nplot=2*setupdataprocessing.nangles;
    else
        nplot=setupdataprocessing.nangles;
    end

    ncol=round((floor(sqrt(nplot*2.5)))/2)*2;
    nrow=ceil(nplot./ncol);
    for iangle=1:length(out)
        if out(iangle).SH.hittingedgeflag
            SHstr='edge';
        else
            SHstr='';
        end

        if setupdataprocessing.nangles==144
            plot_latsum(out(iangle).SH,nrow,ncol,[iangle 0],[num2str(out(iangle).SH.speed,2) ' m/s' SHstr ',' num2str(setupdataprocessing.anglesDeg(iangle)) '^o ' SHstr],[],setupMakePlane.climsforplane(1),setupMakePlane.climsforplane(2),0,0,'w');
        else
            plot_latsum(out(iangle).SH,nrow,ncol,[iangle*2-1 iangle*2],[num2str(out(iangle).SH.speed,2) ' m/s' SHstr],[num2str(iangle) ': ' num2str(setupdataprocessing.anglesDeg(iangle)) '^o ' SHstr],setupMakePlane.climsforplane(1),setupMakePlane.climsforplane(2),0,0,'w');

            for isv=1:length(out(iangle).SecondSW.ilocaltstart)
                plot_latsum(out(iangle).SecondSW,nrow,ncol,[iangle*2-1 iangle*2],[],[],setupMakePlane.climsforplane(1),setupMakePlane.climsforplane(2),1,0,'r',isv);
            end
        end
    end
end



%% reduce saving footprint: don't double-save planes
for i=1:length(out)
    out(i).SH=rmfield(out(i).SH,{'plane','latmm','tms','latsums'});
    out(i).SecondSW=rmfield(out(i).SecondSW,{'latmm','tms','latsums'});
end


end

%% ==========================================================================

function  setupSWS=setRadonSumParams(setupdataprocessing,setupSWS,dataDir)

if regexp(setupdataprocessing.SWSEstimationParams,'_MatchPrevious*')
    % LEGACY: ONE DAY TAKE THIS OUT
    valstring=regexp(setupdataprocessing.SWSEstimationParams,'_MatchPreviouslat(\d+\.?\d?)to(\d+\.?\d?)tmsmin(\d+\.?\d?).*','tokens');
    latstart=str2double(valstring{1}{1});
    latend=str2double(valstring{1}{2});
    tmsminstarttime=str2double(valstring{1}{3});

    for i=1:setupdataprocessing.nangles
        setupSWS.SH(i).minlatmm = latstart;
        setupSWS.SH(i).maxlatmm = latend;
        setupSWS.SH(i).tmsminstarttime=tmsminstarttime; % require wave doesn't start in negative time
        setupSWS.SH(i).maxspeed = 12;
        setupSWS.SH(i).minspeed = 0;
        setupSWS.SH(i).origintmsthreshold=4;
        setupSWS.SH(i).tmsmaxstart=10;
        setupSWS.SH(i).tmsmaxend=20;
    end
    setupSWS.SecondSW=setupSWS.SH;
    return % dont finish running this section, this is overwriting the below
end

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
    setupSWS.SH(i).origintmsthreshold=originthresval;

    % copy global params to where used
    setupSWS.SH(i).minspeed = setupSWS.minspeed;
    setupSWS.SH(i).maxspeed= setupSWS.maxspeed;

    % leave these incase want later, but make super general because somewhat
    % redundant with above
    setupSWS.SH(i).tmsminstarttime=-100; % set to zero to require wave doesn't start in negative time
    setupSWS.SH(i).tmsmaxstart=100;
    setupSWS.SH(i).tmsmaxend=100;

end
setupSWS.SecondSW=setupSWS.SH;


end