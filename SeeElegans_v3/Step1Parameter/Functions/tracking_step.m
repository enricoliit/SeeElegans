function tracking_step(wholestack,outputpath)

% This function tracks neurons with the parameters found in the output
% folder
%
% The inputs are:
%       - wholestack: Stacked recording as a matrix with shape (x,y,z,t)%
%       - outputpath: folder path, specified as a character vector or
%                     string scalar. The folder must contain the files
%                     parameter.mat, parameter2.mat, parameter3.mat and
%                     parameter4.mat

% Example: 
%       >> wholestack = wholestack;
%       >> outputpath = './output/';
%       >> SE_tracking_step(wholestack, outputpath);

% Loading varaibles
load([outputpath '/parameters.mat'],'thresholds','sizes','zsizes','sigmas','voxel_size');
load([outputpath '/parameters3.mat'],'max_gap_closing_distance','max_interframe_distance','non_linking_cost','length_filter');
load([outputpath '/parameters4.mat'],'max_gap_closing_mutual_distance');

% Initialization
addpath('./Functions');
stack_initialization();

maxtime=size(wholestack,4);
number_of_planes=size(wholestack,3);
spotnum_threshold=110;
A_time=[];
max_number_of_tracklets=750;
neurons_tmp=[];
distmatrix=[];
neurons_bv=[];
neurons=[];
neurons_reconstructed_max=[];

% ----------------------------- MAIN PART ---------------------------------
% Step 1: finding spots (LoG filtering and thresholding)
detect_neurons();

% Step 2.1: spot linking through LAP tracking
spot_linking();

% Step 2.2: gap closing based on mutual distance
make_distance_matrix();
track_linking();
validate_tracks();

% Step 3.1: track reconstruction
track_reconstruction();

% Step 3.2: track proximity filter
filter_close_tracks();

% Functions ---------------------------------------------------------------

    % prepare stack for analysis
    function stack_initialization
        wholestack=single(wholestack);
        if size(wholestack,1)>size(wholestack,2)
            wholestack=single(permute(wholestack,[2 1 3 4]));
        end
    end

    % blob detection in each frame and saving
    function detect_neurons
        if ~isfile([outputpath '/A_time.mat'])
            tic
            % finding spots
            disp('Step 1.1: finding spots')
            h = waitbar(0, ['Find spots in ' outputpath]);
            for framenum_id=1:maxtime
                waitbar(framenum_id/maxtime,h);
                spots{framenum_id}=bd3dq(wholestack(:,:,:,framenum_id),thresholds,sizes,zsizes,sigmas,voxel_size,spotnum_threshold);
            end
            elapsed_time_step_1=toc;
            A_time=spots;
            disp(['computation time: ' num2str(elapsed_time_step_1)]);
            save([outputpath '/A_time.mat'],'A_time','elapsed_time_step_1','-v7.3');
            close(h);
        else
            disp('Step 1.1: spots loaded from file');
            load([outputpath '/A_time.mat'],'A_time');
        end
    end

    % LAP tracking
    function spot_linking
        if ~exist('A_time','var')
            load([outputpath '/A_time.mat'],'A_time');
        elseif isempty(A_time)
            load([outputpath '/A_time.mat'],'A_time');
        end
        
        % spot linking ---------------------------------
        if ~isfile([outputpath '/spot_links.mat'])
            disp('Step 1.2: global minimum of interframe spatio-temporal spot proximity');
            spot_links{maxtime}=[];
            disp(['non linking cost = ' num2str(non_linking_cost)]);
            tic
            h = waitbar(0, 'Interframe spot linking ...');
            for t=2:maxtime
                waitbar(t/maxtime,h)
                spot_links{t}=LAP_linker_local_s(t,non_linking_cost,A_time,voxel_size(1),voxel_size(3));
            end
            elapsed_time_step_2=toc;
            save([outputpath '/spot_links.mat'],'spot_links','elapsed_time_step_2','non_linking_cost','-v7.3');
            close(h);
        else
            load([outputpath '/spot_links.mat'],'spot_links');
        end
        
        % building segments lists ---------------------------------
        disp('Step 1.2: building initial segments');
        clear segmentlists segmentlists_time
        [segmentlists_time,segmentlists]=build_links(spot_links);
        
        disp('Step 1.3: gap closing');
        if ~isfile([outputpath '/segment_links.mat'])
            clear segment_links
            segment_links=gap_closer_s(segmentlists_time,segmentlists,max_interframe_distance,max_gap_closing_distance,A_time,voxel_size);
            save([outputpath '/segment_links.mat'],'segment_links','max_interframe_distance','max_gap_closing_distance','-v7.3');
        else
            load([outputpath '/segment_links.mat'],'segment_links');
        end
        
        % building track list
        tracklist=build_tracks(segment_links);
        % filter short linking lists
        disp('trace length filter');
        filtered_exptracelist=merge_and_filter_short_links(tracklist,segmentlists_time,segmentlists,length_filter);
        
        % Creating neurons_tmp variable ---------------------------------
        clear neurons_tmp1 neurons_tmp_unsorted neurons_tmp1_length neurons_tmp
%         neurons_tmp_unsorted(numel(filtered_exptracelist)).Centroid=[];
%         neurons_tmp_unsorted(numel(filtered_exptracelist)).center=[];
        neurons_tmp_unsorted(numel(filtered_exptracelist)).coords=[];
        neurons_tmp_unsorted(numel(filtered_exptracelist)).quality=[];
        
        xy_offset=0;
        z_offset=0;
        h = waitbar(0, 'Building variables...');
        for t=1:maxtime
            waitbar(t/maxtime,h);
            for i=1:numel(filtered_exptracelist)
                if numel(filtered_exptracelist{i})>=t
                    if ~isnan(filtered_exptracelist{i}(t)) && filtered_exptracelist{i}(t)~=0
                        % Centroid_subfield_matrix = cell2mat({A_time{1,t}.Centroid}.');
                        % center_subfield_matrix = cell2mat({A_time{1,t}.center}.');
                        sub_pixel_localization_matrix=A_time{1,t}(:,1:3);
                        sub_pixel_localization_matrix(:,1:2)=sub_pixel_localization_matrix(:,1:2)+xy_offset;
                        sub_pixel_localization_matrix(:,3)=sub_pixel_localization_matrix(:,3)+z_offset;
                        neurons_tmp_unsorted(i).coords=cat(1,neurons_tmp_unsorted(i).coords,single([sub_pixel_localization_matrix(filtered_exptracelist{i}(t),:), t]));
                        % neurons_tmp_unsorted(i).Centroid=cat(1,neurons_tmp_unsorted(i).Centroid,A_time{t}(filtered_exptracelist{i}(t)).Centroid);
                        % neurons_tmp_unsorted(i).center=cat(1,neurons_tmp_unsorted(i).center,A_time{t}(filtered_exptracelist{i}(t),1:3));
                        neurons_tmp_unsorted(i).quality=cat(1,neurons_tmp_unsorted(i).quality,A_time{t}(filtered_exptracelist{i}(t),4));
                    end
                end
            end
        end
        close(h);
        clear filtered_exptracelist
        neurons_tmp=neurons_tmp_unsorted;
        neurons_tmp1_length(numel(neurons_tmp_unsorted))=0;
        for i=1:numel(neurons_tmp_unsorted)
            neurons_tmp1_length(i)=sum(~isnan(neurons_tmp_unsorted(i).coords(:,1)));
        end
        [~,i2]=sort(neurons_tmp1_length,'descend');
        neurons_tmp(numel(neurons_tmp1_length))=[];
        for i=1:numel(neurons_tmp_unsorted)
            neurons_tmp(i)=neurons_tmp_unsorted(i2(i));
        end
        neurons_tmp(max_number_of_tracklets:end)=[];
        save([outputpath '/neurons_tmp.mat'],'neurons_tmp','-v7.3');
    end

    % creating distance matrix
    function make_distance_matrix
        if ~isfile([outputpath '/distmatrix.mat'])
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
            save([outputpath '/distmatrix.mat'],'distmatrix','-v7.3');
        end
    end

    % gap closing based on mutual distance
    function track_linking
        disp('gap closing through mutual distances')
        mutual_position_based_links=[];
        
        if ~exist('distmatrix','var')
            load([outputpath '/distmatrix.mat'],'distmatrix');
        end
        clear links
        distmatrix_mean=squeeze(mean(distmatrix,3,'omitnan'));
        % distmatrix_std=squeeze(std(distmatrix,0,3,'omitnan'));
        taken=[];
        for i=1:numel(neurons_tmp)
            if neurons_tmp(i).coords(end,4) < maxtime
                link=0;
                time_constraint=neurons_tmp(i).coords(end,4);
                loop_set=1:size(distmatrix,1); loop_set(loop_set==i)=[]; loop_set(ismember(loop_set,taken))=[];
                best_candidate=Inf;
                % check=2;
                saved=NaN*ones(1,max(loop_set));
                for j=loop_set
                    if neurons_tmp(j).coords(1,4) >= time_constraint
                        % HO MESSO LA PARENTESI INCLUDENDO NEL QUADRATO LA
                        % VOXEL SIZE
                        current_difference=sqrt( (voxel_size(1).*(distmatrix_mean(i,:,1)-distmatrix_mean(j,:,1))).^2 + (voxel_size(1).*(distmatrix_mean(i,:,2)-distmatrix_mean(j,:,2))).^2 + (voxel_size(3).*(distmatrix_mean(i,:,3)-distmatrix_mean(j,:,3))).^2);
                        current_difference=( mean(current_difference,'omitnan') );
                        saved(j)=current_difference;
                    end
                    [best_candidate, link]=min(saved);
                end
                if link~=0 && best_candidate < max_gap_closing_mutual_distance
                    mutual_position_based_links(i,[1 2 3])=[i,link, best_candidate]; %#ok<AGROW>
                    taken=[taken, link]; %#ok<AGROW>
                end
            end
        end
        skip=[];
        if ~isempty(mutual_position_based_links)
            links_backup=mutual_position_based_links;
            % ex_indices = find(mutual_position_based_links(:,3)>max_gap_closing_mutual_distance);
            % mutual_position_based_links(ex_indices,:) = [];
            mutual_position_based_links(mutual_position_based_links(:,3)>max_gap_closing_mutual_distance,:) = [];
            mutual_position_based_links=mutual_position_based_links(any(mutual_position_based_links,2),:);
            clear reconstruction reconstructed_links
            a=0;
            for i=1:size(mutual_position_based_links,1)
                if mutual_position_based_links(i,1)>0
                    a=a+1;
                    reconstructed_links{a}=[ mutual_position_based_links(i,1) mutual_position_based_links(i,2)]; %#ok<AGROW>
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
                    reconstructed_links{a}=[reconstructed_links{a}, links_backup(found,2)]; %#ok<AGROW>
                end
            end
        end
        clear new_traces
        if exist('reconstructed_links','var')
            new_traces{numel(reconstructed_links)}=[];
            for i=1:numel(reconstructed_links)
                new_traces{i}.coords=[];
%                 new_traces{i}.Centroid=[];
%                 new_traces{i}.center=[];
                new_traces{i}.quality=[];
                for j=1:numel(reconstructed_links{i})
                    new_traces{i}.coords=[new_traces{i}.coords; neurons_tmp(reconstructed_links{i}(j)).coords];
                    % new_traces{i}.Centroid=cat(1,new_traces{i}.Centroid,neurons_tmp(reconstructed_links{i}(j)).Centroid);
%                     new_traces{i}.center=cat(1,new_traces{i}.center,neurons_tmp(reconstructed_links{i}(j)).center);
                    new_traces{i}.quality=cat(1,new_traces{i}.quality,neurons_tmp(reconstructed_links{i}(j)).quality);
                end
            end
            
            % merging linked and complete traces
            skip=unique([reconstructed_links{:}]);
        end
        a=1;
        % clear global neurons_bv
        neurons_bv = neurons_tmp(~ismember(1:numel(neurons_tmp), skip));

        if exist('new_traces','var')
            for i=1:numel(new_traces)
                neurons_bv(a)=new_traces{i};
                a=a+1;
            end
        end
        
        % saving data
        clear distmatrix distmatrix_mean distmatrix_std a best_candidate check ...
            current_difference found i ex_indices j j_id k keep_trace last link ...
            link links links_backup loop_set new_traces reconstructed_links ...
            saved skip taken time_id % framenum time_constraint
        save([outputpath '/neurons_bv.mat'],'neurons_bv','max_gap_closing_mutual_distance');
        
        % visualize
        ncolors=different_colors(numel(neurons_bv));
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
        % if ~exist('neurons_bv','var')
        %     load([outputpath '/neurons_bv.mat'],'neurons_bv');
        % end
    end

    % validation                
    function validate_tracks
        coordinates={neurons_bv.coords};
        lengths=cellfun(@(v)size(v,1),coordinates);
        [~,length_sorting]=sort(lengths,'descend');
        validating=zeros(1,numel(neurons_bv));
        validating(length_sorting)=1;
        save([outputpath '/validation.mat'],'validating');
        a=1;
        clear neurons
        for i=1:numel(validating)
            if validating(i)==1
                neurons(a)=neurons_bv(i); % neurons before validation
                a=a+1;
            end
        end
        clear a
        save([outputpath '/neurons.mat'],'neurons');
    end

    % track reconstruction <------------------- sono arrivato qui a controllare la scala 
    function track_reconstruction
        if ~isfile([outputpath '/neurons_reconstructed_max_1.mat'])
        hf=figure;
        neurons_full(numel(neurons)).coords=[];
        for i=1:numel(neurons)
            neurons_full(i).coords=NaN.*ones(maxtime,4);
            for j=1:size(neurons(i).coords,1)
                if ~isnan(neurons(i).coords(j,4))
                    neurons_full(i).coords(neurons(i).coords(j,4),:)=neurons(i).coords(j,:);
                end
            end
        end
        visualizzati=1:size(neurons,2);
        
        % Step 3.2: track reconstruction
        a=1;
        background_selection=8:36;                                              % SCALA BACKGROUND
        gaussfiltervalue=1;                                                     % SCALA

        looking_window = round(0.85*sizes);                                     % SCALA
        looking_window_small = round(0.6*sizes);                                % SCALA
        radius = sizes;                                                         % SCALA
        radius2= round(0.7*sizes);                                              % SCALA
        radius3= round(0.6*sizes);                                              % SCALA

        hw=waitbar(0,'Measuring traces...');
        for visualizzati_id=1:numel(visualizzati)
            id=visualizzati(visualizzati_id);
            plane=round(mean(neurons(id).coords(:,3)));
            plane(plane<1)=1; plane(plane>number_of_planes)=number_of_planes;
            try
                waitbar(visualizzati_id./numel(visualizzati),hw,['iter = ' num2str(visualizzati_id) '; id = ' num2str(id) '; p = ' num2str(plane)]);
            end
            estimated_mean_position=NaN.*ones(maxtime,3);
            estimated_mean_position(1,:)=neurons_full(id).coords(1,1:3);
            
            % forward
            for time_id=2:maxtime
                if isnan(neurons_full(id).coords(time_id,1))
                    anchor_candidates=[];
                    nearest_candidates=[];
                    for i_id=1:numel(neurons_full)
                        if ~isnan(neurons_full(i_id).coords(time_id,1)) && ~isnan(neurons_full(i_id).coords(time_id-1,1))
                            anchor_candidates=[anchor_candidates, i_id]; %#ok<AGROW>
                            dist_x=voxel_size(1).*((estimated_mean_position(time_id-1,1) - neurons_full(i_id).coords(time_id-1,1)));
                            dist_y=voxel_size(1).*((estimated_mean_position(time_id-1,2) - neurons_full(i_id).coords(time_id-1,2)));
                            dist_z=voxel_size(3).*((estimated_mean_position(time_id-1,3) - neurons_full(i_id).coords(time_id-1,3)));
                            nearest_candidates(i_id,:)=[i_id,sqrt( (dist_x)^2 + (dist_y)^2 + (dist_z)^2 )]; %#ok<AGROW>  % SCALA (differisce fra z e xy, ma forse in questo caso non serve perché parliamo di ricostruzione)
                        end
                    end
                    if ~isempty(nearest_candidates)
                        nearest_candidates=sortrows(nearest_candidates,2);
                    end
                    
                    %                     remove=[];
                    %                     for i=1:size(nearest_candidates,1)
                    %                         if nearest_candidates(i,1)==0
                    %                             remove=[remove, i];
                    %                         end
                    %                     end
                    %                     nearest_candidates(remove,:)=[]; %#ok<AGROW>
                    if ~isempty(nearest_candidates)
                        nearest_candidates((nearest_candidates(:, 1) == 0),:)=[]; %#ok<AGROW>
                        anchor_constraint=20;
                        try %#ok<TRYNC>
                            nearest_candidates([1,anchor_constraint:end],:)=[];
                        end
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
                            dist_x=voxel_size(1).*((estimated_mean_position(time_id+1,1) - neurons_full(i_id).coords(time_id+1,1)));
                            dist_y=voxel_size(1).*((estimated_mean_position(time_id+1,2) - neurons_full(i_id).coords(time_id+1,2)));
                            dist_z=voxel_size(3).*((estimated_mean_position(time_id+1,3) - neurons_full(i_id).coords(time_id+1,3)));
                            nearest_candidates(i_id,:)=[i_id,sqrt( (dist_x)^2 + (dist_y)^2 + (dist_z)^2 )]; % SCALA (differisce fra z e xy, serve perché altrimenti i primi vicini sono tutti i neuroni in quel punto lungo z)
                        end
                    end
                    if ~isempty(nearest_candidates)
                        nearest_candidates=sortrows(nearest_candidates,2);
                    end
                    
                    %                     remove=[];
                    %                     for i=1:size(nearest_candidates,1)
                    %                         if nearest_candidates(i,1)==0
                    %                             remove=[remove, i];
                    %                         end
                    %                     end
                    %                     nearest_candidates(remove,:)=[];
                    
                    if ~isempty(nearest_candidates)
                        nearest_candidates((nearest_candidates(:, 1) == 0),:)=[]; %#ok<AGROW>
                        anchor_constraint=20;
                        try
                            nearest_candidates([1,anchor_constraint:end],:)=[];
                        end
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
            
            % ------------------ SONO ARRIVATO QUI-------------------------
            
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
            
            for i=1:delta_t:(finaltime_traces-3)
                row=round(estimated_mean_position(i,1));
                col=round(estimated_mean_position(i,2));
                plane=round(estimated_mean_position(i,3)); plane(plane<1)=1; plane(plane>number_of_planes)=number_of_planes;
                if ~isnan(row)
                    
                    % ROI AROUND NEURON ---------------------------------------
                    row1=row-looking_window_small; row2=row+looking_window_small; col1=col-looking_window_small; col2=col+looking_window_small;
                    row1(row1<1)=1; row2(row2<1)=1; col1(col1<1)=1; col2(col2<2)=1;
                    row1(row1>s1)=s1; row2(row2>s1)=s1; col1(col1>s2)=s2; col2(col2>s2)=s2;
                    filtered_image=imgaussfilt(wholestack(row1:row2,col1:col2,plane,i),gaussfiltervalue);
                    [row_max,col_max]=find(filtered_image==max(max(filtered_image)));
                    row_1i=row-looking_window; row_2i=row+looking_window; col_1i=col-looking_window; col_2i=col+looking_window;
                    row_1i(row_1i<1)=1; col_1i(col_1i<1)=1;
                    row_2i(row_2i>s1)=s1; col_2i(col_2i>s2)=s2;
                    row_max1=row_1i; row_max2=row_2i;
                    col_max1=col_1i; col_max2=col_2i;
                    row_max1(row_max1<1)=1; row_max2(row_max2<1)=1; col_max1(col_max1<1)=1; col_max2(col_max2<1)=1;
                    row_max1(row_max1>s1)=s1; row_max2(row_max2>s1)=s1; col_max2(col_max2>s2)=s2; col_max1(col_max1>s2)=s2;
                    background=mean(mean(wholestack(background_selection,background_selection,plane,i)));
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
                        ROI_u=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane+1,i:i+3),3)-background)./1,gaussfiltervalue);
                        ROI_u_cP_in=ROI_u(cP_in);
                    catch
                        ROI_u_cP_in=[];
                        ROI_u=[];
                    end
                    ROI_m=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane,i:i+3),3)-background)./1,gaussfiltervalue);
                    try
                        ROI_d=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane-1,i:i+3),3)-background)./1,gaussfiltervalue);
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
            
            % display output
            figure(hf)
            subplot(5,10,a)
            set(gca, 'ColorOrderIndex', 1) %
            trel=numel(trace_reconstruction(1:delta_t:end));
            plot(1:delta_t:trel,(trace_reconstruction(1:delta_t:trel)));
            hold on
            plot(1:delta_t:i,(trace_reconstruction_max(1:delta_t:i)));
            hold on
            plot(1:delta_t:i,(trace_reconstruction_max_no_corona(1:delta_t:i)));
            xlim([0 maxtime]);
            title(num2str(id));
            a=a+1;
            if a>50
                a=1;
                drawnow
                hf=figure;
            end

            % saving variables
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
        end
        close(hw)
        disp('Saving...');
        save([outputpath '/neurons_reconstructed_max_1.mat'],'neurons_reconstructed_max','neurons_reconstructed','-v7.3');
        disp('Saved');
    else
        load([outputpath '/neurons_reconstructed_max_1.mat'],'neurons_reconstructed_max');
    end
    end

    % track proximity filter
    function filter_close_tracks
        disp('check proximity of tracks')
        max_neuron_overlap_xy=sizes/(2*voxel_size(1));  % SCALA
        max_neuron_overlap_z=ceil(zsizes/(2*voxel_size(3)));  % SCALA
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

        % best neurons selection
        traces_to_be_deleted=[];
        quality_check=zeros(1,numel(neurons_reconstructed_max));
        % quality filter
        for i=1:size(too_close_pairs,1)
            for j=1:size(too_close_pairs,2)
                evaluating_neuron=too_close_pairs(i,j);
                already_checked=0;
                if quality_check(evaluating_neuron)~=0
                    already_checked=1;
                end
                if already_checked==0
                    quality_over_time=neurons_reconstructed_max(evaluating_neuron).coords(:,5);
                    % quality averaged over time
                    quality_check(evaluating_neuron)=mean(quality_over_time,'omitnan');
                end
            end
        end
        % removing redundant traces
        a=1;
        for i=1:size(too_close_pairs,1)
            [~, minpos]=min([quality_check(too_close_pairs(i,1)) quality_check(too_close_pairs(i,2))]);
            traces_to_be_deleted(a)=too_close_pairs(i,minpos); %#ok<AGROW>
            a=a+1;
        end
        a=1;
        for i=1:numel(neurons_reconstructed_max)
            if ~ismember(i,traces_to_be_deleted)
                neurons_reconstructed_max_2(a)=neurons_reconstructed_max(i); %#ok<AGROW>
                a=a+1;
            end
        end

        a=1;
        std_threshold=2; % SCALA
        for i=1:numel(neurons_reconstructed_max_2)
            if std(neurons_reconstructed_max_2(i).coords(:,3),0,'omitnan') < std_threshold
                neurons_reconstructed_max_3(a)=neurons_reconstructed_max_2(i); %#ok<AGROW>
                a=a+1;
            end
        end
        clear neurons_reconstructed_max_2
        neurons_reconstructed_max_2=neurons_reconstructed_max_3;
        clear neurons_reconstructed_max_3
        disp('proximity of tracks checked');
        save([outputpath '/neurons_reconstructed_max_2.mat'],'neurons_reconstructed_max_2','thresholds','sizes','sigmas','-v7.3');
        
    end

end


