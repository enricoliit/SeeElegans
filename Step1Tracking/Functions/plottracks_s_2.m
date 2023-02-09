function plottracks_s_2
global display_matrix
global ax
global text_spot_num
global framenum
global plotted_neurons
global colormap_s
a=1;
try
    delete(plotted_neurons);
end
plotted_neurons=[];
display_vector=squeeze(display_matrix(:,framenum,1:2)); % 2D
for i=1:size(display_vector)
    hold on
    plotted_neurons(a)=plot(ax,display_vector(i,2),display_vector(i,1),'o','color',colormap_s(i,:),'MarkerSize',18,'LineWidth',2);
    a=a+1;
end


drawnow
set(text_spot_num,'String',['tracks = ' num2str(a-1)]);

