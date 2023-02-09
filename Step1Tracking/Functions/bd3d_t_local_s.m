function spots=bd3d_t_local_s(wholestack_local,sizes_local,thresholds_local,sigmas_local,zsizes_local)

% volume LOG filtering
filtersize=[sizes_local sizes_local zsizes_local];
h=-fspecial3('log',filtersize,sigmas_local); % this can be taken out of the function
If=imfilter(wholestack_local,h,'symmetric');

% regional maxima detection
regional_max_frame=imregionalmax(If,26);
try
    max_in_A_times_If=regional_max_frame.*wholestack_local;
catch
    nat=whos('wholestack_local');
    class_of_frame=nat.class;
    eval(['regional_max_frame=' class_of_frame '(regional_max_frame);']);
    max_in_A_times_If=regional_max_frame.*wholestack_local; 
end

linear_maxima=find(max_in_A_times_If>thresholds_local);
[x_i,y_i,z_i]=ind2sub(size(wholestack_local),linear_maxima);
qI_i=max_in_A_times_If(linear_maxima);

spots_unfiltered.center=[];
for i=1:numel(x_i)
    spots_unfiltered(i).center=[x_i(i),y_i(i),z_i(i)];
    spots_unfiltered(i).quality=qI_i(i);
end

%% check proximity
overlap_of_neurons=NaN*ones(numel(spots_unfiltered));
a=1;
for i=1:numel(spots_unfiltered)
   for j=i+1:numel(spots_unfiltered)
       mx_i=mean(spots_unfiltered(i).center(1),'omitnan');
       my_i=mean(spots_unfiltered(i).center(2),'omitnan');
       mx_j=mean(spots_unfiltered(j).center(1),'omitnan');
       my_j=mean(spots_unfiltered(j).center(2),'omitnan');
       mz_i=mean(spots_unfiltered(i).center(3),'omitnan');
       mz_j=mean(spots_unfiltered(j).center(3),'omitnan');
       overlap_of_neurons(i,j)=sqrt( (mx_i-mx_j).^2 + (my_i-my_j).^2 );
       if overlap_of_neurons(i,j)<round(sizes_local./2) && abs(mz_i-mz_j)<round(zsizes_local./2)
           too_close_pairs(a,:)=[i, j];
           a=a+1;
       end
   end
end
clear mx_i mx_j my_i my_j

%% remove redundant traces
try
    a=1;
    for i=1:size(too_close_pairs,1)
        [~, minpos]=min([spots_unfiltered(too_close_pairs(i,1)).quality spots_unfiltered(too_close_pairs(i,2)).quality]);
        traces_to_be_deleted(a)=too_close_pairs(i,minpos);
        a=a+1;
    end
end

if ~exist('traces_to_be_deleted','var')
    traces_to_be_deleted=[];
    spots=spots_unfiltered;
else
    a=1;
    for i=1:numel(spots_unfiltered)
        if ~ismember(i,traces_to_be_deleted)
            spots(a)=spots_unfiltered(i);
            a=a+1;
        end
    end
end









