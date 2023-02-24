function gifify(figurehandle,flagfirstframe,filename,delaytime)
% helper script to save animation to .gif file.
% figurehandle is handle of figure to write
% firstflagframe is set to 1 to reset .gif file, any other number will append to current file (optional parameter, default = 0)
% filename is the name of the .gif to be saved (optional parameter, default = testAnimated.gif)
% delaytime is an optional parameter setting speed of .gif (optional parameter, default = 0.1)
%
% use like:
% f1=figure();
% for i=1:100;
%     plot(x(1:i),y(1:i));
%     gifify(f1,i,'test.gif')
% end

if nargin<1
    flagfirstframe=0;
elseif nargin<3
    filename = 'testAnimated.gif';
    delaytime=0.1;
elseif nargin<4
    delaytime=0.1;
end

%axis tight manual % this ensures that getframe() returns a consistent size
drawnow
% Capture the figure as an image
frame = getframe(figurehandle);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if flagfirstframe == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
else % literally anything else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delaytime);
end

