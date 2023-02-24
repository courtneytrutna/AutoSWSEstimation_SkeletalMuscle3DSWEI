function out=CleanAndCategorizeMultiwaves(out,allowflips)
% This function is used after the Radon Sum method to clean up the
% detection of multiple peaks, removing peaks that are not 'independent enough'
% and ultimately outputing at max two SWS estimates (main and one secondary
% wave)

% if allowflips==1, this script will also potentially swap waves
% to assign the 'main wave' as the wave determined to most likely be the SH-wave

% a few goals:
% remove SH wave from the "SV wave" list
% require "SV"'s not to pass/cross SH wave

% note throughout, this script refers to multiwaves as specifically SV
% waves, based on the theory of TI materials, but may infact be handling
% multiple waves from other sources.

out.SH.SHSVflag=''; % add flag to track that swapped

% if only found one wave, assume SH and remove from SV's
if length(out.SecondSW.ilocaltstart)==1
    out.SecondSW.ilocaltstart=[];
    out.SecondSW.ilocaltend=[];
    out.SecondSW.tms_start=[];
    out.SecondSW.tms_end=[];
    out.SecondSW.speed=[];
    return
end

%% remove SV's that are substanitally similar to SH wave

% determine 'corridor' of SH wave
%coded to allow for SV's and SH's to have different lat ranges, but will
%work if they have the same one

tms_start_guess_at_SVlat=out.SH.tms_start+(1/out.SH.speed).*(out.SecondSW.latmm(out.SecondSW.ilatstart)-out.SH.latmm(out.SH.ilatstart));
tms_end_guess_at_SVlat=out.SH.tms_end+(1/out.SH.speed).*(out.SecondSW.latmm(out.SecondSW.ilatend)-out.SH.latmm(out.SH.ilatend));

corridorthresholdstart=1; %ms;
corridorthresholdend=1; %ms
corridorthresholdspeed=.1; % ms
% remove 'SV' waves that fall within corridor
idxremove=[];
for iSV=1:length(out.SecondSW.ilocaltstart)
    test1=abs(out.SecondSW.tms_start(iSV)-tms_start_guess_at_SVlat)<corridorthresholdstart;
    test2=abs(out.SecondSW.tms_end(iSV)-tms_end_guess_at_SVlat)<corridorthresholdend;
    test3=abs(out.SecondSW.speed(iSV)-out.SH.speed)<corridorthresholdspeed;

    if (test1+test2+test3)>=2 % if two out of three are true
        idxremove(length(idxremove)+1)=iSV;
    end
end

if isempty(idxremove) % if no waves fell in coordidor, the acquisition is probalby a mess and the SV's found are noise, so toss ALL of them out
    out.SecondSW.ilocaltstart=[];
    out.SecondSW.ilocaltend=[];
    out.SecondSW.tms_start=[];
    out.SecondSW.tms_end=[];
    out.SecondSW.speed=[];
    out.SecondSW.datapeak=[];
    return

else %otherwise, just remove the SH-like SV wave(s)
    out.SecondSW.ilocaltstart(idxremove)=[];
    out.SecondSW.ilocaltend(idxremove)=[];
    out.SecondSW.tms_start(idxremove)=[];
    out.SecondSW.tms_end(idxremove)=[];
    out.SecondSW.speed(idxremove)=[];
    out.SecondSW.datapeak(idxremove)=[];

end


% and recategorize if SH & SV are switched
if allowflips
    [out]=FlipSHSV(out);
end

% check to see if waves are 'crossing' really early-- usually happens with
% 'merged' wave cases

idxremove=[];
for iSV=1:length(out.SecondSW.ilocaltstart)

    % find equation for SH through all lat
    % y = b*x+y0 = b*(x-x1)+y1
    % tms = 1/speed*(latmm-latstart)+tmsatstart
    latmm_testintercept=2:0.2:20; % if intercept in this range, problem. Note doesn't start at 0, because if intercept in 0-2mm, usually okay

    SH_tmsalllat=out.SH.tms_start+(1/out.SH.speed).*(latmm_testintercept-out.SH.latmm(out.SH.ilatstart));
    SV_tmsalllat=out.SecondSW.tms_start(iSV)+(1/out.SecondSW.speed(iSV)).*(latmm_testintercept-out.SecondSW.latmm(out.SecondSW.ilatstart));

    % look where lines are closest together
    diffSHSVtms=SH_tmsalllat-SV_tmsalllat;
    [~,ilat_closest]=min(abs(diffSHSVtms));

    if ilat_closest==1 % if ilat==1, then don't cross in this lat range, but pre-zero. Assume good.
        continue
    elseif ilat_closest==length(latmm_testintercept) % if closest at furtherst lat range, then SH_speed>SV_speed (true SH and SV, though might be misidentified at this stage). This is bad and we're gonna toss it out
        idxremove(length(idxremove)+1)=iSV;
    else
        % this means they cross at some point, and should be removed
        idxremove(length(idxremove)+1)=iSV;

    end
end

out.SecondSW.ilocaltstart(idxremove)=[];
out.SecondSW.ilocaltend(idxremove)=[];
out.SecondSW.tms_start(idxremove)=[];
out.SecondSW.tms_end(idxremove)=[];
out.SecondSW.speed(idxremove)=[];
out.SecondSW.datapeak(idxremove)=[];

idxremove=[];
%remove SV waves that are still after SH waves
for iSV=1:length(out.SecondSW.ilocaltstart)
    if and(out.SecondSW.ilocaltstart(iSV)>out.SH.itstart,out.SecondSW.ilocaltend(iSV)>out.SH.itend)
        idxremove(length(idxremove)+1)=iSV;
    end
end
out.SecondSW.ilocaltstart(idxremove)=[];
out.SecondSW.ilocaltend(idxremove)=[];
out.SecondSW.tms_start(idxremove)=[];
out.SecondSW.tms_end(idxremove)=[];
out.SecondSW.speed(idxremove)=[];

% reduce to only one SV wave-- first, which will be highest amplitude
if length(out.SecondSW.ilocaltstart)>1
    out.SecondSW.ilocaltstart=out.SecondSW.ilocaltstart(1);
    out.SecondSW.ilocaltend=out.SecondSW.ilocaltend(1);
    out.SecondSW.tms_start=out.SecondSW.tms_start(1);
    out.SecondSW.tms_end=out.SecondSW.tms_end(1);
    out.SecondSW.speed=out.SecondSW.speed(1);
end

% now that removed crossing and picked best, some more might apply to be flipped
if allowflips
    [out]=FlipSHSV(out);
end

end


function [out]=FlipSHSV(out)
if length(out.SecondSW.ilocaltstart)>=1
    if and(out.SecondSW.tms_start(1)>out.SH.tms_start,out.SecondSW.tms_end(1)>out.SH.tms_end) %"SH" is ahead of "SV"
        % when this happens, usually because SV wave was primary detected wave in original
        % SH-lat-range. We then remove the SV wave above, which is
        % incorrect.

        if out.SH.datapeak(1)*.25>out.SecondSW.datapeak(1) % if strongest SV data way smaller than SH data, not actually SV wave, don't switch
            return
        end

        idxwavepull=1; %assume we want to switch the first SV wave
        tmpSH=out.SecondSW;
        tmpSV=out.SH;

        %reassign all relevant values

        % real SH wave is in SV_SHlatrange
        out.SH.avg=[];out.SH.std=[]; % these are based on lateral range, clear out
        out.SH.datapeak  = tmpSH.datapeak(idxwavepull);
        out.SH.speed     = tmpSH.speed(idxwavepull);
        out.SH.itstart   = tmpSH.ilocaltstart(idxwavepull);
        out.SH.itend     = tmpSH.ilocaltend(idxwavepull);
        out.SH.ilatstart = tmpSH.ilatstart;
        out.SH.ilatend   = tmpSH.ilatend;
        out.SH.latsums   = tmpSH.latsums;
        %plane, latmm, tms dont change
        out.SH.tms_end = tmpSH.tms_end(idxwavepull);
        out.SH.tms_start = tmpSH.tms_start(idxwavepull);
        out.SH.hittingedgeflag=0; % clear because it's not for this wave-- shouldn't be here anyway if hit edge
        out.SH.ttpqualmetric=NaN; % clear because it's not for this wave
        if strcmp(out.SH.SHSVflag,'swappedwithSV')
            out.SH.SHSVflag='swappedwithSVx2';
        else
            out.SH.SHSVflag='swappedwithSV'; % add flag to track that swapped
        end

        % real SV wave is the 'SH-like' wave from SV (highest peak)
        out.SecondSW.datapeak     = tmpSV.datapeak;
        out.SecondSW.ilocaltstart = tmpSV.itstart;
        out.SecondSW.ilocaltend   = tmpSV.itend;
        out.SecondSW.ilatstart    = tmpSV.ilatstart;
        out.SecondSW.ilatend      = tmpSV.ilatend;
        out.SecondSW.latsums      = tmpSV.latsums;
        %latmm and tms are the same
        out.SecondSW.tms_end      = tmpSV.tms_end;
        out.SecondSW.tms_start    = tmpSV.tms_start;
        slope = (out.SecondSW.tms_end-out.SecondSW.tms_start)/(out.SecondSW.ilatend-out.SecondSW.ilatstart);
        out.SecondSW.speed = mean(diff(tmpSV.latmm))./slope;
    end
end
end

