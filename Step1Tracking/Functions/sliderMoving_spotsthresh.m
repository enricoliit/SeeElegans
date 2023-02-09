function sliderMoving_spotsthresh
global thresholds
global sizes
global sigmas
global framenum
global slice
global text_threshold
global text_size
global text_sigma
global text_timepoint
global text_slice
global text_spot_num
global I
global wholestack
global ax

set(text_threshold,'String',['threshold = ' num2str(thresholds)]);
set(text_size,'String',['size = ' num2str(sizes)]);
set(text_sigma,'String',['sigma = ' num2str(sigmas)]);
set(text_timepoint,'String',['time = ' num2str(framenum)]);
set(text_slice,'String',['slice = ' num2str(slice)]);
set(text_spot_num,'String',['slice = ' num2str(slice)]);

I=wholestack(:,:,slice,framenum);
ax; imagesc(I); axis equal, axis off;
spots=bd3d;
plotspots(spots)
end