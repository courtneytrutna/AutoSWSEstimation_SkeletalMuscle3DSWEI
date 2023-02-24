function out=plot3Dhelper(angledeg,latmm,dataplot,plottype)
% helper file to plot 3D data in polar-like manner
% angledeg is angles of input data
% latmm is lateral positions of input data
% dataplot is the data to be plotted (size length(angledeg) x
% length(latmm))

% plottype is a flag for different plotting types. The four options are:
% 'points'--> plot dataplot as points in 3D graph, with colors corresponding to lateral position
% 'surf'  --> plot dataplot as surface in 3D graph
% 'k.' or 'grey.' --> plot dataplot as points in 3D graph, with color black (k) or grey dots
% 'grid'  --> will plot polar axes on z=0 plane, giving illusion of polar plotting while using MATLAB's 3D cartesian plotting commands.
% Note 'grid' can be called without inputting data information (ex: plot3Dhelper([],[],[],'grid')

if nargin<4
    plottype='points';
end

if strcmp(plottype,'surf')
    nangle=length(angledeg);
    angledeg(nangle+1)=angledeg(1); % force close surf
    dataplot(nangle+1,:)=dataplot(1,:);
end

[latmesh,anglemesh]=meshgrid(latmm,angledeg);
[x,y,z]=pol2cart(anglemesh/180*pi,latmesh,dataplot);
grid off

switch plottype
    case 'points'
        scatter3(x(:),y(:),z(:),[],latmesh(:),'filled')
    case 'k.'
        out=scatter3(x(:),y(:),z(:),[],'k','filled');
    case 'grey.'
        scatter3(x(:),y(:),z(:),[],[.5 .5 .5],'.')
    case 'surf'
        out=surf(x,y,z,latmesh);
        colormap spring
    case 'grid'
        % plot polar axes on otherwise cartensian plot
        colorline=[.6 .6 .6];
        plot3([-20 20],[0 0],[0 0],'-','Color',colorline)
        plot3([0 0],[-20 20],[0 0],'-','Color',colorline)
        plot3([0 0],[-20 20],[0 0],'-','Color',colorline)

        plot3([-20 20],[-20 20],[0 0],'-','Color',colorline)
        plot3([-20 20],[20 -20],[0 0],'-','Color',colorline)
        plot3([-20 20],[20 -20],[0 0],'-','Color',colorline)

        text(20,0,0,'0^o','Color',colorline)

        for iz=dataplot
            plot3([-20 20],[20 20],[iz iz],'-','Color',colorline)
            plot3([-20 -20],[-20 20],[iz iz],'-','Color',colorline)
        end

        for ir=latmm
            plot3(ir*cos(linspace(0,2*pi,100)),ir*sin(linspace(0,2*pi,100)),0*linspace(0,2*pi,100),'-','Color',colorline)
            text(1,-ir+1,0,[num2str(ir) ' mm'],'Color',colorline)
        end
end
end