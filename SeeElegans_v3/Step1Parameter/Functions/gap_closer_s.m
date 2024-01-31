function links=gap_closer_s(segmentlists_time,segmentlists,max_interframe_distance,max_gap_closing_distance,spots_s,voxel_size)

xy_ratio=voxel_size(1);
z_ratio=voxel_size(3);

% time matrix
segment_start_times=cellfun(@(cell_variable) cell_variable(1), segmentlists_time);
segment_end_times=cellfun(@(cell_variable) cell_variable(end), segmentlists_time);
segment_start_times_label=cellfun(@(cell_variable) cell_variable(1), segmentlists);
segment_end_times_label=cellfun(@(cell_variable) cell_variable(end), segmentlists);
for i=1:numel(segment_start_times)
    for k=1:numel(segment_end_times)
        TC(k,i)=(segment_start_times(i)-segment_end_times(k));
    end
end
TC(TC<=0)=Inf;
TC(TC>max_interframe_distance)=Inf;
TC(TC<Inf)=0;

% Subpixel centroid extraction in previous frame (t1-1), (xc_t0,yc_t0,zc_t0)
for i=1:numel(segment_start_times)
    sub_pixel_localization_matrix_traces_starts(i,:)=spots_s{segment_start_times(i)}(segment_start_times_label(i),1:3);
end
xc_t0=xy_ratio .* sub_pixel_localization_matrix_traces_starts(:,1);
yc_t0=xy_ratio .* sub_pixel_localization_matrix_traces_starts(:,2);
zc_t0=z_ratio .* sub_pixel_localization_matrix_traces_starts(:,3);

for i=1:numel(segment_end_times)
    sub_pixel_localization_matrix_traces_ends(i,:)=spots_s{segment_end_times(i)}(segment_end_times_label(i),1:3);
end
xc_t1=xy_ratio .* sub_pixel_localization_matrix_traces_ends(:,1);
yc_t1=xy_ratio .* sub_pixel_localization_matrix_traces_ends(:,2);
zc_t1=z_ratio .* sub_pixel_localization_matrix_traces_ends(:,3);

% cost matrix calculation
for k = 1:numel(xc_t0)
    for h = 1:length(xc_t1)
        calculated_cost=vecnorm([xc_t0(k), yc_t0(k), zc_t0(k)]-[xc_t1(h), yc_t1(h), zc_t1(h)])+TC(h,k);
        if calculated_cost < max_gap_closing_distance
            C(h,k) = calculated_cost;
        else
            C(h,k) = Inf;
        end
    end
end
C=C.^2;
highest_cost=max(C(C~=Inf));

if isempty(highest_cost)
    highest_cost=100;
end

% lowest cost solution
[links, uR, uC]=matchpairs(single(C),single(highest_cost*1.05));

% add unmatched segments
[gc,grps]=groupcounts([uR; uC]);
unmatched=grps(gc==2);
links=cat(1,links,[unmatched, zeros(numel(unmatched),1)]);