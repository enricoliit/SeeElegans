function spots_filtered=bd3dq(single_volume,thresholds,sizes_local,zsizes_local,sigmas_local,voxel_size,max_spotnum)

sizes=round(sizes_local./voxel_size(1));
zsizes=round(zsizes_local./voxel_size(3));
sigmas=round(sigmas_local./voxel_size(1));
filtersize=[sizes sizes zsizes];

% volume LOG filtering
h=-fspecial3('log',filtersize,sigmas);
If=imfilter(single_volume,h,'symmetric');

% regional maxima detection
regional_max_frame=imregionalmax(If,26);
try
    max_in_single_volume_times_If = regional_max_frame .* single_volume;
catch
    regional_max_frame = cast(regional_max_frame, class(single_volume));
    max_in_single_volume_times_If = regional_max_frame .* single_volume;
end

% quality of spots
linear_maxima = find(max_in_single_volume_times_If > thresholds);
q_i = If(linear_maxima);
qI_i = q_i .* single_volume(linear_maxima);

% quality filter
qorder=sort(qI_i,'descend');
try
    qthresh=qorder(max_spotnum);
    filterout=(qI_i<qthresh);
    linear_maxima(filterout)=[];
    q_i(filterout)=[];
    qI_i(filterout)=[];
end

% spot selection
[x_i,y_i,z_i]=ind2sub(size(single_volume),linear_maxima);
spots_unfiltered = [x_i, y_i, z_i, qI_i];

% proximity check
overlap_of_neurons=NaN*ones(numel(spots_unfiltered(:,1)),numel(spots_unfiltered(:,1)));
a=1;
for i=1:size(spots_unfiltered,1)
   for j=i+1:size(spots_unfiltered,1)
       overlap_of_neurons(i,j)=sqrt( (spots_unfiltered(i,1)-spots_unfiltered(j,1)).^2 + (spots_unfiltered(i,2)-spots_unfiltered(j,2)).^2 );
       if overlap_of_neurons(i,j)<round(sizes./2) && abs(spots_unfiltered(i,3)-spots_unfiltered(j,3))<round(zsizes./2)
           too_close_pairs(a,:)=[i, j];
           a=a+1;
       end
   end
end

% remove redundant tracks
try
    a=1;
    for i=1:size(too_close_pairs,1)
        [~, minpos]=min([qI_i(too_close_pairs(i,1)) qI_i(too_close_pairs(i,2))]);
        traces_to_be_deleted(a)=too_close_pairs(i,minpos);
        a=a+1;
    end
end

if ~exist('traces_to_be_deleted','var')
    spots_filtered = spots_unfiltered;
else
    spots_filtered = spots_unfiltered(~ismember(1:size(spots_unfiltered, 1), traces_to_be_deleted), :);
end