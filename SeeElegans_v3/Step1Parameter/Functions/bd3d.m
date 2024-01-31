function spots=bd3d(single_volume,thresholds,scaled_sizes,scaled_zsizes,scaled_sigmas,voxel_size)

sizes=round(scaled_sizes./voxel_size(1));
zsizes=round(scaled_zsizes./voxel_size(3));
sigmas=round(scaled_sigmas./voxel_size(1));

% volume LOG filtering
h = -fspecial3('log', [sizes sizes zsizes], sigmas);
If = imfilter(single_volume, h, 'symmetric');

% regional maxima detection
regional_max_frame = imregionalmax(If,26);
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

% spot selection
[x_i,y_i,z_i] = ind2sub(size(single_volume), find(max_in_single_volume_times_If > thresholds));
spots_unfiltered = [x_i,y_i,z_i];

% proximity check
% qui si può mettere pdist probabilmente
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
    spots=spots_unfiltered;
else
    spots = spots_unfiltered(~ismember(1:size(spots_unfiltered, 1), traces_to_be_deleted), :);
end