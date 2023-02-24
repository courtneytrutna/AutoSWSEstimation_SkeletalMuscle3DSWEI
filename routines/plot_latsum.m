function axishandle=plot_latsum(out,subplotx,subploty,subplotnum,titlestring1,titlestring2,color_min,color_max,secondplotflag,labelsonflag,colorplot,ipeak)
% helper file to plot output of Radon Sum processing
% out contains the Radon Sum information
% subplotx and subploty are the number of rows and columns respectively

% subplotnum is a two element vector for the subplot number of the
% space-time plane and radon transform space plots. pass in a 0 as the
% second element in subplotnum to not plot the radon transform space

% titlestring1 and titlestring2 are labels for the two subplots
% color_min and color_max are the colorscale used for plotting the space time planes. Radon transform space is allowed to vary
% secondplotflag set to 1 is used to add additional radon trajectories to space-time plane. Use 0 to reset subplot
% labelsonflag set to 1 will use axis labels, set to 0 will not label axes
% colorplot is color of trajectory line and radon transform point
% ipeak is which index of the out structure to plot (can only be single value)


if nargin==8
    secondplotflag = 1;
    labelsonflag=1;
    colorplot='w';
    ipeak=1;
elseif nargin == 9
    labelsonflag=1;
    colorplot='w';
    ipeak=1;
elseif nargin==10
    colorplot='w';
    ipeak=1;
elseif nargin==11
    ipeak=1;
end

%figure(fignum)
%    clf
axishandle(1)=subplot(subplotx,subploty,subplotnum(1));
if ~secondplotflag
    imagesc(out.tms,out.latmm,out.plane)
end
hold on
if ~strcmp(colorplot,'none')
    plot([out.tms_start(ipeak) out.tms_end(ipeak)],[out.latmm(out.ilatstart) out.latmm(out.ilatend)],colorplot, 'LineWidth',1)
end
hold off
if labelsonflag
    xlabel('time (ms)');
    ylabel('lateral position (mm)');
    colorbar
else
    tmp=gca;
    tmp.XTick=[];
    tmp.YTick=[];
end
if ~secondplotflag
    title(titlestring1,'FontSize',8)
end
caxis([color_min color_max])
axis square

if subplotnum(2)
    axishandle(2)=subplot(subplotx,subploty,subplotnum(2));
    if ~secondplotflag
        imagesc(out.tms(1:(end-1)),out.tms(1:(end-1)),out.latsums)
    end
    hold on
    if ~strcmp(colorplot,'none')
        plot([out.tms_end(ipeak) out.tms_end(ipeak)],[out.tms_start(ipeak) out.tms_start(ipeak)],['.' colorplot],'MarkerSize',6)
    end
    hold off
    set(gca,'DataAspectRatio',[1 1 1])
    if ~secondplotflag
        title(titlestring2,'FontSize',8)
    end
    if ~labelsonflag
        tmp=gca;
        tmp.XTick=[];
        tmp.YTick=[];
    end
    %xlabel('end time (ms)')
    %ylabel('start time (ms)')
end