
%==========================================================================

function [out3DSWS,setup3DSWS]=FindSWS_EllipseFit(dataDir,setupdataprocessing,outSWS,setup3DSWS)
% fit an ellipse to the SWS data using a RANSAC method
% dataDir is the location of the data and is used for labeling and saving plots
% setupdataprocessing, outSWS and setup3DSWS are structures with certain expected fields
% setupdataprocessing contains some instructions for fitting
% setup3DSWS contains other instructions for fitting, and is also modified with additional info and passed back out for record keeping
% outSWS contians SWS estimations at each rotation angle
% output variable out3DSWS saves information about ellipse fitting

%% settings for fit

% set outlier threshold
if isempty(setupdataprocessing.setrandsampoutlierthres)
    fitparams.randsampoutlierthreshold=1;
else
    tmp=regexp(setupdataprocessing.setrandsampoutlierthres,'_randsampoutthres(.*)','tokens');
    fitparams.randsampoutlierthreshold=str2double(tmp{1}{1});
end

fitparams.fixedphi = "vary"; % look over 360 for the fiber rotation angle
fitparams.fixedthetatilt=0; % assume fibers are not tilted

%% setup SWS for fitting
% initialize and populate gSWS vector
anglesDeg = setupdataprocessing.anglesDeg;
nangles = length(anglesDeg);
SH_edgeflag   = NaN(1,nangles);
valsqualmetric   = NaN(1,nangles);
gSWS_SHallmeasured   = NaN(1,nangles);
gSWS_SecondWave = NaN(1,nangles);

for iang = 1:nangles
    gSWS_SHallmeasured(iang)   = outSWS(iang).SH.speed;
    SH_edgeflag(iang)=outSWS(iang).SH.hittingedgeflag;
    if isempty(outSWS(iang).SH.qualmetric)
        valsqualmetric(iang)=0;
    else
        valsqualmetric(iang)=outSWS(iang).SH.qualmetric;
    end
    gSWS_SecondWave(1:length(outSWS(iang).SecondSW.speed),iang) = outSWS(iang).SecondSW.speed; % sometimes has multiple, will resize accordingly
    gSWS_SecondWave(max(length(outSWS(iang).SecondSW.speed),2):size(gSWS_SecondWave,1),iang)=NaN; % fill with NaN instead of zero. Sometimes length=0, so then need max(length,2). first row already zeros from initiatlization, only need for 2+
    valsqualmetric_SecondWave(1:length(outSWS(iang).SecondSW.speed),iang)=0; %don't carry this information through
end


SH_edgeflag(isnan(SH_edgeflag))=0;
SH_edgeflag=logical(SH_edgeflag);

%remove bad data
gSWS_SH=gSWS_SHallmeasured;
gSWS_SH(SH_edgeflag)   = NaN; %remove bad values

%% Preform Fit and Organize Output

%run RANSAC ellipse fitting
[outFitSWSvsAngleSH,removedoutlierflag]=FitSWSvsAngle_RANSAC(gSWS_SH,anglesDeg,fitparams.fixedphi,fitparams.fixedthetatilt,fitparams,setup3DSWS);

outFitSWSvsAngleSH.removedoutlierflag=removedoutlierflag; % Save for future plotting

outFitSWSvsAngleSH.percentptsincluded=100*sum(~isnan(outFitSWSvsAngleSH.cvalsToFit))./length(outFitSWSvsAngleSH.cvalsToFit);
outFitSWSvsAngleSH.gSWS_SHallmeasured=gSWS_SHallmeasured;
outFitSWSvsAngleSH.gSWS_SVallmeasured=gSWS_SecondWave;

% Final outvalues: consistent format across all methods (cPar,cPerp) plus
% holding on to information of interest

out3DSWS.cPar=outFitSWSvsAngleSH.cParfit;
out3DSWS.cPerp=outFitSWSvsAngleSH.cPerpfit;
out3DSWS.phiRot=outFitSWSvsAngleSH.phifit;

out3DSWS.percentpts=outFitSWSvsAngleSH.percentptsincluded;
out3DSWS.costfunctionval=outFitSWSvsAngleSH.ellipsecostfunctionval;
out3DSWS.costfunctionval_inlier=outFitSWSvsAngleSH.ellipsecost_inlier;

out3DSWS.EachSWSEstimate=outSWS;
out3DSWS.RANSACEllipse=outFitSWSvsAngleSH;
setup3DSWS.RANSACellipsefitparams=fitparams;

%% Plotting

if setup3DSWS.fignum
    % pull relevant parameters for plotting
    anglesfullSH = outFitSWSvsAngleSH.anglesfull;
    cfitfullSH   = outFitSWSvsAngleSH.cfitfull;

    tmp=regexp(dataDir,'/','split');titlestring=tmp{end};titlestring=replace(titlestring,'_',' ');

    if sum(isnan(valsqualmetric))==length(valsqualmetric)
        valsqualmetric(:)=0; % overwrite nans with zeros so below works
    end

    % Add ellipse fit to plots from each sub-acquistion

    tmp=figure(1);
    tmp=findobj(tmp,'Type','Axes');
    if length(tmp)>2
        % move everything over a bit
        for i=1:length(tmp);tmp(i).Position(1)=tmp(i).Position(1)-.1;end

        % mark excluded
        try % this might not always work depending on SWS method
            for i=find(removedoutlierflag)
                angletmp=setupdataprocessing.anglesDeg(i);
                for i2=1:length(tmp)
                    findaxes(i2)=strcmp([num2str(i) ': ' num2str(angletmp) '^o '],tmp(i2).Title.String);
                end
                iaxmodify=find(findaxes);
                tmp(iaxmodify).Title.String=[tmp(iaxmodify).Title.String 'X'];
            end
        catch
            disp('Warning: exclusions not marked in subplot titles')
        end

        % add polar plot in top left
        polax=polaraxes;
        polax.Position=[.78 .65 .25 .25];
        if ~sum(valsqualmetric,'omitnan')==0 %if not all zeros, plot quality metric under outlier
            polarscatter(anglesDeg(outFitSWSvsAngleSH.removedoutlierflag).*pi/180,gSWS_SHallmeasured(outFitSWSvsAngleSH.removedoutlierflag),20,valsqualmetric(removedoutlierflag),'filled');hold on;
            for i=1:size(gSWS_SecondWave,1); polarscatter(anglesDeg*pi/180,gSWS_SecondWave(i,:)',20,valsqualmetric_SecondWave(i,:),'^','filled');end
            for i=1:size(gSWS_SecondWave,1); polarscatter(anglesDeg*pi/180,gSWS_SecondWave(i,:)',20,'r','^');end
        else
            for i=1:size(gSWS_SecondWave,1); polarscatter(anglesDeg*pi/180,gSWS_SecondWave(i,:)',5,'r','filled');hold on;end

        end
        polarplot(anglesDeg(removedoutlierflag).*pi/180,gSWS_SHallmeasured(removedoutlierflag),'kx') % plot removed outliers
        hold on
        polarplot(anglesDeg(SH_edgeflag).*pi/180,gSWS_SHallmeasured(SH_edgeflag),'rx') % plot radons excluded due to edgehitting
        polarscatter(outFitSWSvsAngleSH.anglesToFit*pi/180,outFitSWSvsAngleSH.cvalsToFit,20,valsqualmetric,'filled')
        caxis([0 1])
        rlim([0 8])
        polarplot(anglesfullSH*pi/180,cfitfullSH,'k-') %plot SWS data from fit

        %plot fiber angle
        polarplot([outFitSWSvsAngleSH.phifit.*pi/180,outFitSWSvsAngleSH.phifit*pi/180],[0 6],'k--')

        title({['inlier MSE:' num2str(outFitSWSvsAngleSH.ellipsecost_inlier)],...
            ['included: ' num2str(outFitSWSvsAngleSH.percentptsincluded,3) '%'],...
            ['total ellipse cost function val: ' num2str(outFitSWSvsAngleSH.ellipsecostfunctionval) ]})

        %plot SWS values
        ax_vals=axes;
        ax_vals.Position=[.85 .4 .1 .15];
        offset1=0;
        offset2=0.1;
        plot((1:2)-offset1-offset2,[outFitSWSvsAngleSH.cParfit,outFitSWSvsAngleSH.cPerpfit],'k*')
        legend('SH fit','Location','southoutside')
        tmp=gca;tmp.XTick=[1,2];tmp.XTickLabel={'c_{Par}','c_{Perp}'}; %label axis
        xlim([.5 2.5]);ylim([0 10])
        ylabel('m/s')
        title({['SWS_L= ' num2str(outFitSWSvsAngleSH.cParfit,3) '; SWS_T= ' num2str(outFitSWSvsAngleSH.cPerpfit,3)]})
        ax_vals.Position=[.85 .4 .1 .15]; %reset after legend

        % plot angles
        ax_deg=axes;
        ax_deg.Position=[.85 .1 .1 .15];
        plot(1-offset1-offset2,outFitSWSvsAngleSH.phifit,'k*')
        tmp=gca;tmp.XTick=[];
        xlim([.5 1.5]);ylim([-90 90])
        ylabel('Rot Angle (degrees)')
        legendstrings={'SH fit'};


        legend(legendstrings,'Location','southoutside')
        if isstruct(fitparams.fixedthetatilt)
            title({['Rot =  ' num2str(outFitSWSvsAngleSH.phifit) '^o'], ['TiltEllipse = ' num2str(fitparams.fixedthetatilt.tiltforellipse) '^o, TiltcL = ' num2str(fitparams.fixedthetatilt.tiltcorrectionforCL) '^o']})
        else
            title(['Rot =  ' num2str(outFitSWSvsAngleSH.phifit) '^o, Tilt = ' num2str(fitparams.fixedthetatilt) '^o'])
        end
        ax_deg.Position=[.85 .1 .1 .15]; %reset after legend

        sgtitle(titlestring)

        [savefolder,SWSsettingsname]=GenerateSaveFileName(dataDir,setupdataprocessing);
        print('-dpng',[savefolder 'SpaceTimeAndEllipsePlots' SWSsettingsname '.png'])
    end


end

end
