function SE_tracking_step(wholestack,outputpath)
% thi function tracks neurons with the parameters found in the outputfolder
% The inputs are:
%       - wholestack: Stacked recording as a matrix with shape (x,y,z,t)%
%       - outputpath: folder path, specified as a character vector or
%                     string scalar. The folder must contain the files
%                     parameter.mat, parameter2.mat, parameter3.mat and
%                     parameter4.mat
%
% Example: SE_tracking_step(wholestack, outputpath)
%       - wholestack = wholestack;
%       - outputpath = './output';

%% Initialization
addpath('./Functions');
% initial parameters
global xy_pixel_size; global z_pixel_size
% parameters to be chosen
global thresholds; global sizes; global zsizes; global sigmas; global time_processing
% visualization
global framenum; global I; global slice; global ax
global wholestack; global video_length
global number_of_planes
% on screen text
global text_threshold; global text_size; global text_sigma; global text_timepoint;
global text_slice; global text_spot_num;
% results
global A_time; global maxtime; global neurons_bv_full;
global non_linking_cost
global max_gap_closing_distance
global max_interframe_distance
global length_filter
global max_gap_closing_mutual_distance
xy_pixel_size=1;
z_pixel_size=2.*xy_pixel_size;
% initial varaibles
thresholds=200; sizes=12; zsizes=3; sigmas=3; framenum=500; slice=1; time_processing=1;

if ~exist(outputpath, 'dir')
    mkdir(outputpath)
end

load([outputpath '\parameters.mat']);
load([outputpath '\parameters2.mat']);
load([outputpath '\parameters3.mat']);
load([outputpath '\parameters4.mat']);
if 1==1
    wholestack=single(wholestack);
    for i=1:12
        wholestack(:,:,i,:)=Kalman_Stack_Filter(squeeze(wholestack(:,:,i,:)),0.7,1);
    end
    if size(wholestack,1)>size(wholestack,2)
        wholestack=single(permute(wholestack,[2 1 3 4]));
    end
end
number_of_planes=size(wholestack,3);
video_length=size(wholestack,4);
maxtime=video_length;

%% Step 1: finding spots (LoG filtering and thresholding)
local_maxtime=maxtime;
sizes_local=sizes;
thresholds_local=thresholds;
sigmas_local=sigmas;
zsizes_local=zsizes;
spotnum_threshold=110;

if ~isfile([outputpath '\A_time.mat'])
    tic
    % finding spots
    disp('Step 1.1: finding spots')
    spots{local_maxtime}=[];
    h = waitbar(0, ['Find spots in ' stackpath]);
    for framenum_id=1:local_maxtime
        waitbar(framenum_id/local_maxtime,h);
        spots{framenum_id}=bd3d_t_local_filtered(wholestack(:,:,:,framenum_id),sizes_local,thresholds_local,sigmas_local,zsizes_local,spotnum_threshold);
    end
    elapsed_time_step_1=toc;
    A_time=spots;
    disp(['computation time: ' num2str(elapsed_time_step_1)]);
    save([outputpath '\A_time.mat'],'A_time','elapsed_time_step_1','-v7.3');
    close(h);
else
    disp('Step 1.1: spots loaded from file');
    load([outputpath '\A_time.mat']);
end
%% LAP tracking: spot linking

if ~exist('A_time','var')
    load([outputpath '\A_time.mat']);
elseif isempty(A_time)
    load([outputpath '\A_time.mat']);
end
xy_pixel_size_local=xy_pixel_size;
z_pixel_size_local=z_pixel_size;
if ~isfile([outputpath '\spot_links.mat'])
    clear spot_links
    disp('Step 1.2: global minimum of interframe spatio-temporal spot proximity');
    spot_links{maxtime}=[];
    disp(['non linking cost = ' num2str(non_linking_cost)]);
    tic
    h = waitbar(0, 'Interframe spot linking ...');
    for t=2:maxtime
        waitbar(t/maxtime)
        spot_links{t}=LAP_linker_local(t,non_linking_cost,A_time,xy_pixel_size_local,z_pixel_size_local);
    end
    elapsed_time_step_2=toc;
    save([outputpath '\spot_links.mat'],'spot_links','elapsed_time_step_2','non_linking_cost','-v7.3');
    close(h);
else
    load([outputpath '\spot_links.mat'],'spot_links');
end
% building segments lists
disp('Step 1.3: building initial segments');
clear segmentlists segmentlists_time
[segmentlists_time,segmentlists]=build_links(spot_links);

%% LAP tracking: gap closing

disp('Step 1.3: gap closing');
if ~isfile([outputpath '\segment_links.mat'])
    clear segment_links
    segment_links=gap_closer(A_time,segmentlists_time,segmentlists,max_interframe_distance,max_gap_closing_distance,xy_pixel_size_local,z_pixel_size_local);
    save([outputpath '\segment_links.mat'],'segment_links','max_interframe_distance','max_gap_closing_distance','-v7.3');
else
    load([outputpath '\segment_links.mat'],'segment_links');
end
% building track list
tracklist=build_tracks(segment_links);
% filter short linking lists
disp('trace length filter');
filtered_exptracelist=merge_and_filter_short_links(tracklist,segmentlists_time,segmentlists,length_filter);

%% Creating neurons_tmp variable

clear neurons_tmp1 neurons_tmp_unsorted neurons_tmp1_length neurons_tmp
neurons_tmp_unsorted(numel(filtered_exptracelist)).Centroid=[];
neurons_tmp_unsorted(numel(filtered_exptracelist)).center=[];
neurons_tmp_unsorted(numel(filtered_exptracelist)).quality=[];
neurons_tmp_unsorted(numel(filtered_exptracelist)).coords=[];
xy_offset=0;
z_offset=0;
for t=1:maxtime
    for i=1:numel(filtered_exptracelist)
        if ~isnan(filtered_exptracelist{i}(t)) && filtered_exptracelist{i}(t)~=0
            %
            Centroid_subfield_matrix = cell2mat({A_time{1,t}.Centroid}.');
            center_subfield_matrix = cell2mat({A_time{1,t}.center}.');
            sub_pixel_localization_matrix=center_subfield_matrix;
            sub_pixel_localization_matrix(:,1:2)=sub_pixel_localization_matrix(:,1:2)+xy_offset;
            sub_pixel_localization_matrix(:,3)=sub_pixel_localization_matrix(:,3)+z_offset;
            neurons_tmp_unsorted(i).coords=cat(1,neurons_tmp_unsorted(i).coords,single([sub_pixel_localization_matrix(filtered_exptracelist{i}(t),:), t]));
            neurons_tmp_unsorted(i).Centroid=cat(1,neurons_tmp_unsorted(i).Centroid,A_time{t}(filtered_exptracelist{i}(t)).Centroid);
            neurons_tmp_unsorted(i).center=cat(1,neurons_tmp_unsorted(i).center,A_time{t}(filtered_exptracelist{i}(t)).center);
            neurons_tmp_unsorted(i).quality=cat(1,neurons_tmp_unsorted(i).quality,A_time{t}(filtered_exptracelist{i}(t)).quality);
        end
    end
end
clear filtered_exptracelist
neurons_tmp=neurons_tmp_unsorted;
for i=1:numel(neurons_tmp_unsorted)
    neurons_tmp1_length(i)=sum(~isnan(neurons_tmp_unsorted(i).coords(:,1)));
end
[i1,i2]=sort(neurons_tmp1_length,'descend');
neurons_tmp(numel(neurons_tmp1_length))=[];
for i=1:numel(neurons_tmp_unsorted)
    neurons_tmp(i)=neurons_tmp_unsorted(i2(i));
end
% ATTENZIONE INTERROMPO NEURONS_TMP
neurons_tmp(750:end)=[];
save([outputpath '\neurons_tmp.mat'],'neurons_tmp','-v7.3');

%% Step 2.1: creating distance matrix

disp('creating distance matrix')
clear distmatrix
distmatrix=NaN*ones(numel(neurons_tmp),numel(neurons_tmp),1,3);
for i=1:numel(neurons_tmp)
    valid_times=neurons_tmp(i).coords(:,4);
    logic_times=~isnan(valid_times);
    valid_times_ids=valid_times(logic_times);
    for k=1:numel(neurons_tmp)
        valid_times_k=neurons_tmp(k).coords(:,4);
        logic_times_k=~isnan(valid_times_k);
        valid_times_ids_k=valid_times_k(logic_times_k);
        i_matrix=NaN.*ones(maxtime,3); i_matrix(valid_times_ids,:)=neurons_tmp(i).coords(logic_times,1:3);
        k_matrix=NaN.*ones(maxtime,3); k_matrix(valid_times_ids_k,:)=neurons_tmp(k).coords(logic_times_k,1:3);
        distmatrix(i,k,1,1:3)= mean((i_matrix - k_matrix),1,'omitnan');
    end
end
save([outputpath '\distmatrix.mat'],'distmatrix','-v7.3');

%% Step 2.2: gap closing based on mutual distance
% parameters
disp('gap closing through mutual distances')
mutual_position_based_links=[];

if ~exist('distmatrix','var')
    load([outputpath '\distmatrix.mat'],'distmatrix');
end
clear links
distmatrix_mean=squeeze(mean(distmatrix,3,'omitnan'));
distmatrix_std=squeeze(std(distmatrix,0,3,'omitnan'));
taken=[];
for i=1:numel(neurons_tmp)
    if neurons_tmp(i).coords(end,4)<maxtime
        link=0;
        time_constraint=neurons_tmp(i).coords(end,4);
        loop_set=1:size(distmatrix,1); loop_set(loop_set==i)=[]; loop_set(ismember(loop_set,taken))=[];
        best_candidate=Inf;
        check=2;
        saved=NaN*ones(1,max(loop_set));
        for j=loop_set
            if neurons_tmp(j).coords(1,4)>=time_constraint
                current_difference=sqrt( xy_pixel_size.*(distmatrix_mean(i,:,1)-distmatrix_mean(j,:,1)).^2 + xy_pixel_size.*(distmatrix_mean(i,:,2)-distmatrix_mean(j,:,2)).^2 + z_pixel_size.*(distmatrix_mean(i,:,3)-distmatrix_mean(j,:,3)).^2); % aggiungo il tempo?
                current_difference=( mean(current_difference,'omitnan') );
                saved(j)=current_difference;
            end
            [best_candidate, link]=min(saved);
        end
        if link~=0 && best_candidate<max_gap_closing_mutual_distance
            mutual_position_based_links(i,[1 2 3])=[i,link, best_candidate];
            taken=[taken, link];
        end
    end
end
skip=[];
if ~isempty(mutual_position_based_links)
    links_backup=mutual_position_based_links;
    ex_indices = find(mutual_position_based_links(:,3)>max_gap_closing_mutual_distance);
    mutual_position_based_links(ex_indices,:) = [];
    mutual_position_based_links=mutual_position_based_links(any(mutual_position_based_links,2),:);
    clear reconstruction reconstructed_links
    a=0;
    for i=1:size(mutual_position_based_links,1)
        if mutual_position_based_links(i,1)>0
            a=a+1;
            reconstructed_links{a}=[ mutual_position_based_links(i,1) mutual_position_based_links(i,2)];
            last=mutual_position_based_links(i,2);
            mutual_position_based_links(i,:)=[0 0 0];
            if last==0
                continue
            end
            while last~=0
                found=find(mutual_position_based_links(:,1)==last);
                if numel(found)>0
                    reconstructed_links{a}=[reconstructed_links{a}, mutual_position_based_links(found,2)];
                    last=mutual_position_based_links(found,2);
                    mutual_position_based_links(found,:)=[0 0 0];
                else
                    break
                end
            end
            reconstructed_links{a}=[reconstructed_links{a}, links_backup(found,2)];
        end
    end
end
clear new_traces
if exist('reconstructed_links','var')
    new_traces{numel(reconstructed_links)}=[];
    for i=1:numel(reconstructed_links)
        new_traces{i}.coords=[];
        new_traces{i}.Centroid=[];
        new_traces{i}.center=[];
        new_traces{i}.quality=[];
        for j=1:numel(reconstructed_links{i})
            new_traces{i}.coords=[new_traces{i}.coords; neurons_tmp(reconstructed_links{i}(j)).coords];
            new_traces{i}.Centroid=cat(1,new_traces{i}.Centroid,neurons_tmp(reconstructed_links{i}(j)).Centroid);
            new_traces{i}.center=cat(1,new_traces{i}.center,neurons_tmp(reconstructed_links{i}(j)).center);
            new_traces{i}.quality=cat(1,new_traces{i}.quality,neurons_tmp(reconstructed_links{i}(j)).quality);
        end
    end
    % merging linked and complete traces
    skip=unique([reconstructed_links{:}]);
end
a=1;
clear global neurons_bv
for i=1:size(neurons_tmp,2)
    if ismember(i,skip)==0
        neurons_bv(a)=neurons_tmp(i);
        a=a+1;
    end
end
if exist('new_traces','var')
    for i=1:numel(new_traces)
        neurons_bv(a)=new_traces{i};
        a=a+1;
    end
end


%% saving data
clear distmatrix distmatrix_mean distmatrix_std a best_candidate check ...
    current_difference found i ex_indices j j_id k keep_trace last link ...
    link links links_backup loop_set new_traces reconstructed_links ...
    saved skip taken time_id % framenum time_constraint
save([outputpath '\neurons_bv.mat'],'neurons_bv','max_gap_closing_mutual_distance');
%% visualize
ncolors=distinguishable_colors(numel(neurons_bv));
figure
for i=1:numel(neurons_bv)
    plot(neurons_bv(i).coords(:,2),neurons_bv(i).coords(:,1),'Color',ncolors(i,:))
    hold on
    text(double(neurons_bv(i).coords(1,2)),double(neurons_bv(i).coords(1,1)),num2str(i))
    hold on
end
title('Reconstruction')
set(gca,'YDir','reverse')
axis equal
axis off
clear ncolors i
if ~exist('neurons_bv','var')
    load([outputpath '\neurons_bv.mat']);
end

%% visualize traces
lower=[1 2 7 8];
center=[3 4 9 10];
upper=[5 6 11 12];
validating=-1.*ones(1,numel(neurons_bv));
maximum_width=size(wholestack,1);
maximum_height=size(wholestack,2);
coordinates={neurons_bv.center};
lengths=cellfun(@(v)size(v,1),coordinates);
[length_values,length_sorting]=sort(lengths,'descend');
validating=zeros(1,numel(neurons_bv));
validating(length_sorting)=1;
clear lower upper center
save([outputpath '\validation.mat'],'validating');

%% validated neurons
a=1;
clear neurons
for i=1:numel(validating)
    if validating(i)==1
        neurons(a)=neurons_bv(i); % neurons before validation
        a=a+1;
    end
end
clear a
save([outputpath '\neurons.mat'],'neurons');

%% Step 3.1: second distance matrix
clear distmatrix2
distmatrix2=NaN*ones(numel(neurons),numel(neurons),1,3);
timestep=100;
hw=waitbar(0,'distance matrix calculation');
for i=1:numel(neurons)
    waitbar(i/numel(neurons),hw,'distance matrix calculation')
    j_id=1;
    for j=1:timestep:size(neurons(i).coords,1)
        time=neurons(i).coords(j,4);
        for k=1:numel(neurons)
            clear time_id
            time_id=find(neurons(k).coords(:,4)==time);
            try
                time_id=time_id(1);
            end
            if isempty(time_id)==0
                distmatrix2(i,k,j_id,1:3)= (neurons(i).coords(j,1:3) - neurons(k).coords(time_id,1:3)) ;
            end
        end
        j_id=j_id+1;
    end
end
close(hw);
distmatrix2(distmatrix2==0)=NaN;
distmatrix2_mean=squeeze(mean(distmatrix2,3,'omitnan'));
distmatrix2_std=squeeze(std(distmatrix2,0,3,'omitnan'));
save([outputpath '\distmatrix2.mat'],'distmatrix2','-v7.3');
clear i j j_id k t framenum time_id traces
figure;
a=1;
b=1;
visualizzati=[];
for i=1:numel(neurons)
    neurons_full(i).coords=NaN.*ones(maxtime,4);
    for j=1:size(neurons(i).coords,1)
        if ~isnan(neurons(i).coords(j,4))
            neurons_full(i).coords(neurons(i).coords(j,4),:)=neurons(i).coords(j,:);
        end
    end
    visualizzati=[visualizzati, i];
end
visualizzati=1:size(neurons,2);

%% Step 3.2: track reconstruction
a=1;
looking_window=round(sizes/2);
looking_window_small=round(sizes/3);
radius = round(2*sizes/3);
radius2= round(sizes/3+sizes/10);
radius3= round(sizes/3);
w=-fspecial('log',[sizes sizes],sigmas);
if ~exist('distmatrix2_mean','var')
    distmatrix2_mean=(mean(distmatrix2,3,'omitnan'));
    distmatrix_mean=(mean(distmatrix,3,'omitnan'));
end
hw=waitbar(0,'Measuring traces...');
for visualizzati_id=1:numel(visualizzati)
    id=visualizzati(visualizzati_id);
    plane=round(mean(neurons(id).coords(:,3)));
    plane(plane<1)=1; plane(plane>number_of_planes)=number_of_planes;
    waitbar(visualizzati_id./numel(visualizzati),hw,['iter = ' num2str(visualizzati_id) '; id = ' num2str(id) '; p = ' num2str(plane)]);
    estimated_mean_position=NaN.*ones(maxtime,3); % PERCHE' NON METTERLO A NAN ????
    estimated_mean_position(1,:)=neurons_full(id).coords(1,1:3);
    % forward
    for time_id=2:maxtime
        if isnan(neurons_full(id).coords(time_id,1))
            anchor_candidates=[];
            nearest_candidates=[];
            for i_id=1:numel(neurons_full)
                if ~isnan(neurons_full(i_id).coords(time_id,1)) && ~isnan(neurons_full(i_id).coords(time_id-1,1))
                    anchor_candidates=[anchor_candidates, i_id];
                    dist_x=(estimated_mean_position(time_id-1,1) - neurons_full(i_id).coords(time_id-1,1));
                    dist_y=(estimated_mean_position(time_id-1,2) - neurons_full(i_id).coords(time_id-1,2));
                    dist_z=(estimated_mean_position(time_id-1,3) - neurons_full(i_id).coords(time_id-1,3));
                    nearest_candidates(i_id,:)=[i_id,sqrt( (dist_x)^2 + (dist_y)^2 + (dist_z)^2 )];
                end
            end
            if ~isempty(nearest_candidates)
                nearest_candidates=sortrows(nearest_candidates,2);
            end
            remove=[];
            for i=1:size(nearest_candidates,1)
                if nearest_candidates(i,1)==0
                    remove=[remove, i];
                end
            end
            nearest_candidates(remove,:)=[];
            anchor_constraint=20;
            try
                nearest_candidates([1,anchor_constraint:end],:)=[];
            end
            if ~isempty(nearest_candidates)
                anchor_candidates=nearest_candidates(:,1);
            else
                anchor_candidates=[];
            end
            % -------------------------------------------------------------
            estimated_position=0;
            number_of_anchors=0;
            for i_id=1:numel(anchor_candidates)
                i=anchor_candidates(i_id);
                anchor_timepoint_before=time_id-1;
                anchor_timepoint=time_id;
                if isempty(anchor_timepoint)
                    continue
                end
                anchor_timepoint=anchor_timepoint(1);
                dist_x=(estimated_mean_position(anchor_timepoint_before,1) - neurons_full(i).coords(anchor_timepoint_before,1));
                dist_y=(estimated_mean_position(anchor_timepoint_before,2) - neurons_full(i).coords(anchor_timepoint_before,2));
                dist_z=(estimated_mean_position(anchor_timepoint_before,3) - neurons_full(i).coords(anchor_timepoint_before,3));
                if ~isnan(neurons_full(i).coords(anchor_timepoint,1)) && ~isnan(dist_x) && ~isnan(dist_y) && ~isnan(dist_z)
                    addendum=[dist_x, dist_y, dist_z]+neurons_full(i).coords(anchor_timepoint,1:3);
                    if ~isnan(addendum)
                        estimated_position=estimated_position+addendum;
                        number_of_anchors=number_of_anchors+1;
                    end
                end
            end
            estimated_mean_position(time_id,:)=estimated_position./number_of_anchors;
        else
            estimated_mean_position(time_id,:)=neurons_full(id).coords(time_id,1:3);
        end
    end
    % backwards
    for time_id=maxtime-1:-1:1
        if isnan(estimated_mean_position(time_id,1))
            % -------------------------------------------------------------
            anchor_candidates=[];
            nearest_candidates=[];
            for i_id=1:numel(neurons_full)
                if ~isnan(neurons_full(i_id).coords(time_id,1)) && ~isnan(neurons_full(i_id).coords(time_id+1,1))
                    anchor_candidates=[anchor_candidates, i_id];
                    dist_x=(estimated_mean_position(time_id+1,1) - neurons_full(i_id).coords(time_id+1,1));
                    dist_y=(estimated_mean_position(time_id+1,2) - neurons_full(i_id).coords(time_id+1,2));
                    dist_z=(estimated_mean_position(time_id+1,3) - neurons_full(i_id).coords(time_id+1,3));
                    nearest_candidates(i_id,:)=[i_id,sqrt( (dist_x)^2 + (dist_y)^2 + (dist_z)^2 )];
                end
            end
            if ~isempty(nearest_candidates)
                nearest_candidates=sortrows(nearest_candidates,2);
            end
            remove=[];
            for i=1:size(nearest_candidates,1)
                if nearest_candidates(i,1)==0
                    remove=[remove, i];
                end
            end
            nearest_candidates(remove,:)=[];
            anchor_constraint=20;
            try
                nearest_candidates([1,anchor_constraint:end],:)=[];
            end
            if ~isempty(nearest_candidates)
                anchor_candidates=nearest_candidates(:,1);
            else
                anchor_candidates=[];
            end
            % -------------------------------------------------------------
            estimated_position=0;
            number_of_anchors=0;
            for i_id=1:numel(anchor_candidates)
                i=anchor_candidates(i_id);
                anchor_timepoint_before=time_id+1;
                anchor_timepoint=time_id;
                if isempty(anchor_timepoint)
                    continue
                end
                anchor_timepoint=anchor_timepoint(1);
                dist_x=(estimated_mean_position(anchor_timepoint_before,1) - neurons_full(i).coords(anchor_timepoint_before,1));
                dist_y=(estimated_mean_position(anchor_timepoint_before,2) - neurons_full(i).coords(anchor_timepoint_before,2));
                dist_z=(estimated_mean_position(anchor_timepoint_before,3) - neurons_full(i).coords(anchor_timepoint_before,3));
                if ~isnan(neurons_full(i).coords(anchor_timepoint,1)) && ~isnan(dist_x) && ~isnan(dist_y) && ~isnan(dist_z)
                    addendum=[dist_x, dist_y, dist_z]+neurons_full(i).coords(anchor_timepoint,1:3);
                    if ~isnan(addendum)
                        estimated_position=estimated_position+addendum;
                        number_of_anchors=number_of_anchors+1;
                    end
                end
            end
            estimated_mean_position(time_id,:)=estimated_position./number_of_anchors;
        end
    end
    [s1, s2]=size(wholestack,[1 2]);
    clear trace_reconstruction trace_reconstruction_max
    trace_reconstruction=NaN.*ones(1,maxtime);
    trace_reconstruction_max=NaN.*ones(1,maxtime);
    delta_t=1;
    finaltime_traces=maxtime;
    clear trace_ROI selected_pixels
    clear c_row c_col row col row1 row2 col1 col2 row_max col_max row_1i row_2i col_1i col_2i row_max1 row_max2 col_max1 col_max2
    neurons_reconstructed_max(visualizzati_id).coords=NaN.*ones(maxtime,5);
    neurons_reconstructed(visualizzati_id).coords=NaN.*ones(maxtime,5);
    for i=1:delta_t:(finaltime_traces-3)%(numel(neurons_full(id).coords(:,1))-3) % 96 1096
        row=round(estimated_mean_position(i,1));
        col=round(estimated_mean_position(i,2));
        plane=round(estimated_mean_position(i,3)); plane(plane<1)=1; plane(plane>number_of_planes)=number_of_planes;
        if ~isnan(row)
            % ROI AROUND NEURON ---------------------------------------
            row1=row-looking_window_small; row2=row+looking_window_small; col1=col-looking_window_small; col2=col+looking_window_small;
            row1(row1<1)=1; row2(row2<1)=1; col1(col1<1)=1; col2(col2<2)=1;
            row1(row1>s1)=s1; row2(row2>s1)=s1; col1(col1>s2)=s2; col2(col2>s2)=s2;
            filtered_image=imgaussfilt(wholestack(row1:row2,col1:col2,plane,i),1);
            [row_max,col_max]=find(filtered_image==max(max(filtered_image)));
            row_1i=row-looking_window; row_2i=row+looking_window; col_1i=col-looking_window; col_2i=col+looking_window;
            row_1i(row_1i<1)=1; col_1i(col_1i<1)=1;
            row_2i(row_2i>s1)=s1; col_2i(col_2i>s2)=s2;
            row_max1=row_1i; row_max2=row_2i;
            col_max1=col_1i; col_max2=col_2i;
            row_max1(row_max1<1)=1; row_max2(row_max2<1)=1; col_max1(col_max1<1)=1; col_max2(col_max2<1)=1;
            row_max1(row_max1>s1)=s1; row_max2(row_max2>s1)=s1; col_max2(col_max2>s2)=s2; col_max1(col_max1>s2)=s2;
            background=mean(mean(wholestack(8:36,8:36,plane,i)'));
            trace_reconstruction(i)=mean(mean((wholestack(row1:row2,col1:col2,plane,i)-background)./1));
            row_c=row1+row_max-1; col_c=col1+col_max-1;
            [columnsInImage, rowsInImage] = meshgrid(1:s2, 1:s1);
            center_col = col_c;
            center_row = row_c;
            center_col=center_col(1);
            center_row=center_row(1);
            c_col(i,1)=center_col;
            c_row(i,1)=center_row;
            circlePixels = (rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 <= radius.^2 & (rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 >= radius2.^2;
            circlePixels_in = ((rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 <= radius3.^2);
            cP=(circlePixels(row_1i:row_2i,col_1i:col_2i));
            cP_in=(circlePixels_in(row_1i:row_2i,col_1i:col_2i));
            try
                ROI_u=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane+1,i:i+3),3)-background)./1,1);
                ROI_u_cP_in=ROI_u(cP_in);
            catch
                ROI_u_cP_in=[];
                ROI_u=[];
            end
            ROI_m=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane,i:i+3),3)-background)./1,1);
            try
                ROI_d=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane-1,i:i+3),3)-background)./1,1);
                ROI_d_cP_in=ROI_d(cP_in);
            catch
                ROI_d_cP_in=[];
                ROI_d=[];
            end
            pixels_from_ROIs=[ROI_u_cP_in; ROI_m(cP_in); ROI_d_cP_in];
            [chosen_pixels, chosen_pixels_location]=sort(pixels_from_ROIs,'descend');
            threshold_on_pixels=round(numel(chosen_pixels)*0.1)+1;
            chosen_pixels(threshold_on_pixels:end)=[];
            chosen_pixels_location(threshold_on_pixels:end)=[];
            trace_ROI{i}=cat(3,ROI_d,ROI_m,ROI_u);
            selected_pixels{i}=chosen_pixels_location;
            if ~isempty((chosen_pixels))
                if ~isempty(ROI_m(find(cP)))
                    trace_reconstruction_max(i)=median(chosen_pixels)-median(ROI_m(find(cP)));
                else
                    trace_reconstruction_max(i)=NaN;
                end
                trace_reconstruction_max_no_corona(i)=median(chosen_pixels);
            else
                trace_reconstruction_max_no_corona(i)=NaN;
            end
        else
            trace_reconstruction_max(i)=NaN;
            trace_reconstruction_max_no_corona(i)=NaN;
        end
    end
    subplot(5,10,a)
    set(gca, 'ColorOrderIndex', 1) %
    trel=numel(trace_reconstruction(1:delta_t:end));
    p1=plot(1:delta_t:trel,(trace_reconstruction(1:delta_t:trel)));
    hold on
    p12=plot(1:delta_t:i,(trace_reconstruction_max(1:delta_t:i)));
    hold on
    p15=plot(1:delta_t:i,(trace_reconstruction_max_no_corona(1:delta_t:i)));
    xlim([0 maxtime]);
    title(num2str(id));
    valid_times_in_trace_reconstruction=(~isnan(trace_reconstruction));
    valid_times_in_trace_reconstruction_max=(~isnan(trace_reconstruction_max));
    timepoints=1:maxtime; timepoints_max=1:maxtime;
    timepoints=timepoints(valid_times_in_trace_reconstruction);
    timepoints_max=timepoints_max(valid_times_in_trace_reconstruction_max);
    neurons_reconstructed(visualizzati_id).coords(timepoints,5)=trace_reconstruction(timepoints);
    neurons_reconstructed_max(visualizzati_id).coords(timepoints_max,5)=(trace_reconstruction_max(timepoints_max));
    neurons_reconstructed_max(visualizzati_id).coords(timepoints_max,[1 2])=[c_row(timepoints_max,:), c_col(timepoints_max,:)];
    neurons_reconstructed_max(visualizzati_id).coords(timepoints_max,[3 4])=[estimated_mean_position(timepoints_max,3), (timepoints_max)'];
    neurons_reconstructed_max(visualizzati_id).selected_pixels=selected_pixels(cellfun(@(v)~isempty(v),selected_pixels));
    a=a+1;
    if a>50
        a=1;
        drawnow
        figure
    end
end
close(hw)
disp('Saving...');
save([outputpath '\neurons_reconstructed_max_1.mat'],'neurons_reconstructed_max','neurons_reconstructed','-v7.3');
disp('Saved');
%% checking neuron proximity
clear neurons_reconstructed_max_1
load([outputpath '\neurons_reconstructed_max_1.mat'],'neurons_reconstructed_max','neurons_reconstructed');
%%
disp('check proximity of tracks')
max_neuron_overlap_xy=sizes/2;
max_neuron_overlap_z=ceil(zsizes/2);
overlap_of_neurons=NaN*ones(numel(neurons_reconstructed_max));
a=1;
clear too_close_pairs
for i=1:numel(neurons_reconstructed_max)
    for j=i+1:numel(neurons_reconstructed_max)
        mx_i=mean(neurons_reconstructed_max(i).coords(:,1),'omitnan');
        my_i=mean(neurons_reconstructed_max(i).coords(:,2),'omitnan');
        mx_j=mean(neurons_reconstructed_max(j).coords(:,1),'omitnan');
        my_j=mean(neurons_reconstructed_max(j).coords(:,2),'omitnan');
        mz_i=mean(neurons_reconstructed_max(i).coords(:,3),'omitnan');
        mz_j=mean(neurons_reconstructed_max(j).coords(:,3),'omitnan');
        overlap_of_neurons(i,j)=sqrt( (mx_i-mx_j).^2 + (my_i-my_j).^2 );
        if overlap_of_neurons(i,j)<max_neuron_overlap_xy && abs(mz_i-mz_j)<max_neuron_overlap_z
            too_close_pairs(a,:)=[i, j];
            a=a+1;
        end
    end
end
clear mx_i mx_j my_i my_j mz_i mz_j

%% best neurons selection
traces_to_be_deleted=[];
quality_check=zeros(1,numel(neurons_reconstructed_max));
for collapsing=1
    % quality filter
    w=-fspecial('log',[sizes sizes],sigmas);
    for i=1:size(too_close_pairs,1)
        for j=1:size(too_close_pairs,2)
            evaluating_neuron=too_close_pairs(i,j);
            already_checked=0;
            if quality_check(evaluating_neuron)~=0
                already_checked=1;
            end
            if already_checked==0
                quality_over_time=neurons_reconstructed_max(evaluating_neuron).coords(timepoints_max,5);
                % qualità media nel tempo
                quality_check(evaluating_neuron)=mean(quality_over_time,'omitnan');
            end
        end
    end
    % removing redundant traces
    a=1;
    for i=1:size(too_close_pairs,1)
        [~, minpos]=min([quality_check(too_close_pairs(i,1)) quality_check(too_close_pairs(i,2))]);
        traces_to_be_deleted(a)=too_close_pairs(i,minpos);
        a=a+1;
    end
    a=1;
    for i=1:numel(neurons_reconstructed_max)
        if ~ismember(i,traces_to_be_deleted)
            neurons_reconstructed_max_2(a)=neurons_reconstructed_max(i);
            a=a+1;
        end
    end
end

%%
a=1;
std_threshold=2;
for i=1:numel(neurons_reconstructed_max_2)
    std(neurons_reconstructed_max_2(i).coords(:,3),0,'omitnan')
    if std(neurons_reconstructed_max_2(i).coords(:,3),0,'omitnan')<std_threshold
        neurons_reconstructed_max_3(a)=neurons_reconstructed_max_2(i);
        a=a+1;
    end
end
clear neurons_reconstructed_max_2
neurons_reconstructed_max_2=neurons_reconstructed_max_3;
clear neurons_reconstructed_max_3
disp('proximity of tracks checked');
save([outputpath '\neurons_reconstructed_max_2.mat'],'neurons_reconstructed_max_2','thresholds','sizes','sigmas','-v7.3');