function sliderMoving_tracks_s_2
global framenum
global slice
global text_timepoint
global text_slice
global I
global wholestack_s
global shown_fig

set(text_timepoint,'String',['time = ' num2str(framenum)])
set(text_slice,'String',['slice = ' num2str(slice)])

I=wholestack_s(:,:,slice,framenum);
shown_fig.CData=I;
plottracks_s
end