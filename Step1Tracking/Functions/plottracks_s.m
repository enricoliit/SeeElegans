function plottracks_s
global ax
global text_spot_num
global framenum
global A_time_s
global plotted_neurons
global neurons_tmp_unsorted_s
a=1;
try
   delete(plotted_neurons); 
end
plotted_neurons=[];
for i=1:numel(A_time_s{framenum})
    if ~isempty(A_time_s{framenum}(i).color)
        hold on
        plotted_neurons(a)=plot(ax,A_time_s{framenum}(i).center(2),A_time_s{framenum}(i).center(1),'o','color',A_time_s{framenum}(i).color,'MarkerSize',18,'LineWidth',2);
        a=a+1;
    end
end
drawnow
set(text_spot_num,'String',['tracks = ' num2str(a-1) ' / ' num2str(numel(neurons_tmp_unsorted_s)) ]);
end
