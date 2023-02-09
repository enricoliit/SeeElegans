function SE_parameter_selection(wholestack, outputpath)
% thi function opens a GUI to choose tracking parameters to track neurons
% in 4D calcium imaging recordings. 
% The inputs are:
%       - wholestack: Stacked recording as a matrix with shape (x,y,z,t)
%
%       - outputpath: folder path, specified as a character vector or 
%                     string scalar. If the folder does not exist, it will
%                     be created
%
% Example: SE_parameter_selection(wholestack, outputpath)
%       - wholestack = wholestack;
%       - outputpath = './output';


%% Initialization
addpath('./Functions/');

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

% selected parameters
global maxtime;
global non_linking_cost
global max_gap_closing_distance
global max_interframe_distance
global length_filter
global max_gap_closing_mutual_distance

% initial definitions
xy_pixel_size=1;
z_pixel_size=2.*xy_pixel_size;

% initial varaibles
thresholds=200; sizes=12; zsizes=3; sigmas=3; framenum=500; slice=1; time_processing=1;

disp(outputpath)
if ~exist(outputpath, 'dir')
    mkdir(outputpath)
end
number_of_planes=size(wholestack,3);
video_length=size(wholestack,4);
maxtime=video_length;
audio=0;

%% Step 0.0: visualization and segmentation parameter selection
disp('Starting...')
global go;
f=figure('Units','normalized','Position',[0.125 0.125 0.75 0.75]);
set(f, 'MenuBar', 'none');
ax=axes(f,'Position',[0.05 0.2 0.9 0.8]);
framenum(framenum>video_length)=video_length;
I=wholestack(:,:,slice,framenum)-mean(mean(wholestack(1:38,1:38,slice,framenum)));
imagesc(I), axis equal, axis off
% slider range definition
min_threshold=2; max_threshold=max(wholestack(:)); thresholdrange=max_threshold-min_threshold;
minsize=2; maxsize=25; sizerange=maxsize-minsize;
minsigma=1; maxsigma=150; sigmarange=maxsigma-minsigma;
mintime=1; maxtime=video_length; timerange=maxtime-mintime;
minslice=1; maxslice=number_of_planes; slicerange=maxslice-minslice; slicerange(slicerange==0)=1;


guiapp.sld_threshold = uicontrol('parent',f,'Style','slider','SliderStep',[1 25]./1000,'Value',thresholds,'min',min_threshold,'max',max_threshold,'Position',[50 10 150 25]);
guiapp.sld_size = uicontrol('parent',f,'Style','slider','Value',sizes,'min',minsize,'max',maxsize,'SliderStep',[1 5]./sizerange,'Position',[250 10 150 25]);
guiapp.sld_sigma = uicontrol('parent',f,'Style','slider','Value',sigmas,'min',minsigma,'max',maxsigma,'SliderStep',[1 5]./sigmarange,'Position',[450 10 150 25]);
guiapp.sld_timepoint= uicontrol('parent',f,'Style','slider','Value',framenum,'min',mintime,'max',maxtime,'SliderStep',[1 5]./timerange,'Units','normalized','Position',[0.05 0.1525 0.9 0.025]);
guiapp.sld_slice= uicontrol('parent',f,'Style','slider','Value',slice,'min',minslice,'max',maxslice,'SliderStep',[1 5]./slicerange,'Units','normalized','Position',[0.05 0.1275 0.9 0.025]);
guiapp.sld_time_processing = uicontrol('parent',f,'Style','checkbox','String','framenum processing','Value',1,'Position',[850 10 150 25]);
text_threshold = uicontrol('parent',f,'Style','text','string',['threshold = ' num2str(thresholds)],'Position',[50 35 150 25]);
text_size = uicontrol('parent',f,'Style','text','string',['size = ' num2str(sizes)],'Position',[250 35 150 25]);
text_sigma = uicontrol('parent',f,'Style','text','string',['sigma = ' num2str(sigmas)],'Position',[450 35 150 25]);
text_timepoint = uicontrol('parent',f,'Style','text','string',['framenum = ' num2str(framenum)],'Position',[650 35 150 25]);
text_slice = uicontrol('parent',f,'Style','text','string',['slice = ' num2str(slice)],'Position',[650 10 150 25]);
text_spot_num = uicontrol('parent',f,'Style','text','string',['spot number = ' num2str(0)],'Position',[650 60 150 25]);

set(guiapp.sld_threshold,'Callback','global thresholds; thresholds=guiapp.sld_threshold.Value; sliderMoving_spotsthresh;');
set(guiapp.sld_size,'Callback','global sizes; sizes=guiapp.sld_size.Value; sliderMoving_spotsthresh;');
set(guiapp.sld_sigma,'Callback','global sigmas; sigmas=guiapp.sld_sigma.Value; sliderMoving_spotsthresh;');
set(guiapp.sld_time_processing,'Callback', 'global time_processing; time_processing=time_processing-1; sliderMoving_spotsthresh;');
set(guiapp.sld_timepoint,'Callback','global framenum; framenum=round(guiapp.sld_timepoint.Value); sliderMoving_spotsthresh;');
set(guiapp.sld_slice,'Callback','global slice; slice=guiapp.sld_slice.Value; sliderMoving_spotsthresh;');
guiapp.bttn_go= uicontrol('parent',f,'Style','pushbutton','String','GO!','Value',0,'Units','normalized','Position',[0.95 0.0475 0.04 0.05],'Callback','global go; go=1;');
go=0;
colormap('gray')
drawnow

if ~isfile([outputpath '\parameters.mat'])
    while go==0
        drawnow
    end
    save([outputpath '\parameters.mat'],'thresholds','sizes','sigmas','zsizes','time_processing','-v7.3');
else
    load([outputpath '\parameters.mat'])
end
structfun(@delete,guiapp);

%% Step 0.1: visualization and subvideo cropping

global parameters
global texts
global go;
global check;
parameters.first_frame=1;
parameters.last_frame=maxtime;
guiapp.first_frame = uicontrol('parent',f,'Style','pushbutton','Value',1,'Position',[1050 60 125 25],'String','Start Frame');
guiapp.last_frame = uicontrol('parent',f,'Style','pushbutton','Value',550,'Position',[1050 30 125 25],'String','End Frame');
texts.first_frame = uicontrol('parent',f,'Style','text','string',['first frame = ' num2str(parameters.first_frame)],'Position',[1175 60 175 25]);
texts.last_frame = uicontrol('parent',f,'Style','text','string',['last frame = ' num2str(parameters.last_frame)],'Position',[1175 30 175 25]);
set(guiapp.first_frame,'Callback','global parameters; global texts; global framenum; guiapp.first_frame.Value=framenum; parameters.first_frame=guiapp.first_frame.Value; set(texts.first_frame,"string",["first frame = " num2str(parameters.first_frame)]);');
set(guiapp.last_frame,'Callback','global parameters; global texts; global framenum; guiapp.last_frame.Value=framenum; parameters.last_frame=guiapp.last_frame.Value; set(texts.last_frame,"string",["first frame = " num2str(parameters.last_frame)]);');
guiapp.sld_timepoint= uicontrol('parent',f,'Style','slider','Value',framenum,'min',mintime,'max',maxtime,'SliderStep',[1 5]./timerange,'Units','normalized','Position',[0.05 0.1525 0.9 0.025]);
guiapp.sld_slice= uicontrol('parent',f,'Style','slider','Value',slice,'min',minslice,'max',maxslice,'SliderStep',[1 5]./slicerange,'Units','normalized','Position',[0.05 0.1275 0.9 0.025]);
set(guiapp.sld_timepoint,'Callback','global framenum; framenum=round(guiapp.sld_timepoint.Value); sliderMoving_spotsthresh;');
set(guiapp.sld_slice,'Callback','global slice; slice=guiapp.sld_slice.Value; sliderMoving_spotsthresh;');
guiapp.bttn_go= uicontrol('parent',f,'Style','pushbutton','String','GO!','Value',0,'Units','normalized','Position',[0.95 0.0475 0.04 0.05],'Callback','global go; go=1;');
go=0;
global wholestack_s
if ~isfile([outputpath '\parameters2.mat'])
    while go==0
        drawnow
    end
    disp(['Evaluating from ' num2str(parameters.first_frame) ' to ' num2str(parameters.last_frame) '...']);
    wholestack_s=wholestack(:,:,:,parameters.first_frame:parameters.last_frame);
    save([outputpath '\parameters2.mat'],'parameters','-v7.3');
else
    load([outputpath '\parameters2.mat']);
end
wholestack_s=wholestack(:,:,:,parameters.first_frame:parameters.last_frame);
structfun(@delete,guiapp);
structfun(@delete,texts);

%% Step 1: finding spots (LoG filtering and thresholding) in subvideo
global NUMERO_PIANI_s
global DURATA_TOTALE_s
global maxtime_s
global wholestack_s
global A_time_s
global xy_pixel_size; global z_pixel_size
global thresholds; global sizes; global zsizes; global sigmas;
NUMERO_PIANI_s=size(wholestack_s,3);
DURATA_TOTALE_s=size(wholestack_s,4);
maxtime_s=DURATA_TOTALE_s;
local_maxtime=maxtime_s;
sizes_local=sizes;
thresholds_local=thresholds;
sigmas_local=sigmas;
zsizes_local=zsizes;
if ~isfile([outputpath '\A_time_s.mat']) || ~exist('A_time_s','var')
    tic
    % finding spots -----------
    disp('Step 1.1: finding spots')
    spots_s{local_maxtime}=[];
    h = waitbar(0, 'Finding spots ...');
    for framenum_id=1:local_maxtime
        spots_s{framenum_id}=bd3d_t_local_s(wholestack_s(:,:,:,framenum_id),sizes_local,thresholds_local,sigmas_local,zsizes_local);
        waitbar(framenum_id/local_maxtime,h);
    end
    close(h)
    elapsed_time_step_1=toc;
    A_time_s=spots_s;
    disp(['computation time: ' num2str(elapsed_time_step_1)]);
    save([outputpath '\A_time_s.mat'],'A_time_s');
else
    if ~exist('A_time_s','var')
        load([outputpath '\A_time_s.mat'],'A_time_s');
    end
end

%% Choosing tracking parameters to be checked in subvideo
global shown_fig
global ax
clear global plotted_neurons
try
    figure(f)
catch
    f=figure('Units','normalized','Position',[0.125 0.125 0.75 0.75]);
end
% set(f, 'MenuBar', 'none');
ax=axes(f,'Position',[0.05 0.2 0.9 0.8]);
framenum(framenum>DURATA_TOTALE_s)=DURATA_TOTALE_s;
I=wholestack_s(:,:,slice,framenum)-mean(mean(wholestack_s(1:38,1:38,slice,framenum)));
shown_fig=imagesc(I); axis equal, axis off
check=0;
go=0;
non_linking_cost=5;
max_gap_closing_distance=5;
max_interframe_distance=2;
length_filter=1;
dx_texts=200;
if ~isfile([outputpath '\parameters3.mat'])
    while check==0 && go==0
        % LAP tracking: spot linking
        clear spot_links_s filtered_exptracelist_s
        A_time_local=A_time_s; xy_pixel_size_local=xy_pixel_size; z_pixel_size_local=z_pixel_size; spot_links_s{maxtime_s}=[];
        disp('Step 1.2: global minimum of interframe spatio-temporal spot proximity');
        h = waitbar(0, 'Interframe spot linking ...');
        tic
        % LAP linker
        for t=2:maxtime_s
            spot_links_s{t}=LAP_linker_local_s(t,non_linking_cost,A_time_local,xy_pixel_size_local,z_pixel_size_local);
            waitbar(t/maxtime_s,h);
        end
        close(h)
        elapsed_time_step_2=toc;
        % building segments lists
        disp('Step 1.3: building initial segments');
        clear segmentlists_s segmentlists_time_s
        [segmentlists_time_s,segmentlists_s]=build_links(spot_links_s);
        % LAP tracking: gap closing
        disp('Step 1.4: gap closing');
        clear segment_links
        segment_links_s=gap_closer_s(segmentlists_time_s,segmentlists_s,max_interframe_distance,max_gap_closing_distance,A_time_s);
        % building track list
        tracklist_s=build_tracks(segment_links_s);
        % filter short linking lists
        disp('trace length filter');
        filtered_exptracelist_s=merge_and_filter_short_links_s(tracklist_s,segmentlists_time_s,segmentlists_s,length_filter);
        %
        clear global neurons_tmp_unsorted_s colormap_s
        global colormap_s
        global neurons_tmp_unsorted_s
        try
            neurons_tmp_unsorted_s(numel(filtered_exptracelist_s)).coords=[];
        end
        xy_offset=-floor(sizes/2)+1;
        z_offset=-floor(zsizes/2)+1;
        colormap_s=distinguishable_colors(numel(filtered_exptracelist_s));
        for i=1:numel(filtered_exptracelist_s)
            neurons_tmp_unsorted_s(i).coords=NaN.*ones(maxtime_s,4);
        end
        try
            delete(test_s);
        end
        % -------------------------------------------------------------------------
        try
            structfun(@delete,guiapp);
            delete(text_non_linking_cost);
            delete(text_max_gap_closing_distance);
            delete(text_max_interframe_distance);
            delete(text_length_filter);
        end
        check=1;
        shown_fig.CData=I;
        % slider range definition
        maxsize=25; minsize=2; sizerange=maxsize-minsize;
        minsigma=3; maxsigma=150; sigmarange=maxsigma-minsigma;
        mintime=1; maxtime_s=DURATA_TOTALE_s; timerange=maxtime_s-mintime;
        minslice=1; maxslice=number_of_planes; slicerange=maxslice-minslice; slicerange(slicerange==0)=1;
        disp('Tracking parameters selection...')
        guiapp.sld_timepoint= uicontrol('parent',f,'Style','slider','Value',framenum,'min',mintime,'max',maxtime_s,'SliderStep',[1 5]./timerange,'Units','normalized','Position',[0.05 0.1525 0.9 0.025]);
        guiapp.sld_slice= uicontrol('parent',f,'Style','slider','Value',slice,'min',minslice,'max',maxslice,'SliderStep',[1 5]./slicerange,'Units','normalized','Position',[0.05 0.1275 0.9 0.025]);
        text_timepoint = uicontrol('parent',f,'Style','text','string',['framenum = ' num2str(framenum)],'Position',[650+dx_texts 35 150 25]);
        text_slice = uicontrol('parent',f,'Style','text','string',['slice = ' num2str(slice)],'Position',[650+dx_texts 10 150 25]);
        text_spot_num = uicontrol('parent',f,'Style','text','string',['spot number = ' num2str(0)],'Position',[650+dx_texts 60 150 25]);
        set(guiapp.sld_timepoint,'Callback','global framenum; framenum=round(guiapp.sld_timepoint.Value); sliderMoving_tracks_s_2;');
        set(guiapp.sld_slice,'Callback','global slice; slice=guiapp.sld_slice.Value; sliderMoving_tracks_s_2;');
        guiapp.non_linking_cost = uicontrol('parent',f,'Style','edit','String',num2str(non_linking_cost),'Position',[50 20 150 25]);
        guiapp.max_gap_closing_distance = uicontrol('parent',f,'Style','edit','String',num2str(max_gap_closing_distance),'Position',[250 20 150 25]);
        guiapp.max_interframe_distance = uicontrol('parent',f,'Style','edit','String',num2str(max_interframe_distance),'Position',[450 20 150 25]);
        guiapp.length_filter = uicontrol('parent',f,'Style','edit','String',num2str(length_filter),'Position',[650 20 150 25]);
        text_non_linking_cost = uicontrol('parent',f,'Style','text','string',['non linking cost = ' num2str(non_linking_cost)],'Position',[50 45 150 25]);
        text_max_gap_closing_distance = uicontrol('parent',f,'Style','text','string',['max gap closing distance = ' num2str(max_gap_closing_distance)],'Position',[250 45 150 25]);
        text_max_interframe_distance = uicontrol('parent',f,'Style','text','string',['max interframe distance = ' num2str(max_interframe_distance)],'Position',[450 45 150 25]);
        text_length_filter = uicontrol('parent',f,'Style','text','string',['length filter = ' num2str(length_filter)],'Position',[650 45 150 25]);
        set(guiapp.non_linking_cost,'Callback','global non_linking_cost; non_linking_cost=str2num(guiapp.non_linking_cost.String); set(text_non_linking_cost,"string",["non linking cost = " (guiapp.non_linking_cost.String)]);');
        set(guiapp.max_gap_closing_distance,'Callback','global max_gap_closing_distance; max_gap_closing_distance=str2num(guiapp.max_gap_closing_distance.String); set(text_max_gap_closing_distance,"string",["max gap closing distance = " (guiapp.max_gap_closing_distance.String)]);');
        set(guiapp.max_interframe_distance,'Callback','global max_interframe_distance; max_interframe_distance=str2num(guiapp.max_interframe_distance.String); set(text_max_interframe_distance,"string",["max interframe distance = " (guiapp.max_interframe_distance.String)]);');
        set(guiapp.length_filter,'Callback','global length_filter; length_filter=str2num(guiapp.length_filter.String); set(text_length_filter,"string",["length filter = " (guiapp.length_filter.String)]);');
        guiapp.bttn_go= uicontrol('parent',f,'Style','pushbutton','String','GO!','Value',0,'Units','normalized','Position',[0.95 0.0475 0.04 0.05],'Callback','go=1;');
        guiapp.bttn_check= uicontrol('parent',f,'Style','pushbutton','String','Check','Value',0,'Units','normalized','Position',[0.85 0.0475 0.04 0.05],'Callback','check=0;');
        go=0;
        colormap('gray')
        % -------------------------------------------------------------------------
        for i=1:numel(filtered_exptracelist_s)
            for t=1:maxtime_s
                if ~isnan(filtered_exptracelist_s{i}(t)) && filtered_exptracelist_s{i}(t)~=0
                    %
                    id=filtered_exptracelist_s{i}(t);
                    sub_pixel_localization=cell2mat({A_time_s{1,t}(id).center}.');
                    sub_pixel_localization=sub_pixel_localization;
                    neurons_tmp_unsorted_s(i).coords(t,:)=single([sub_pixel_localization, t]);
                    id=filtered_exptracelist_s{i}(t);
                    A_time_s{1,t}(id).color=colormap_s(i,:);
                end
            end
            hold on
            test_s(i)=plot(neurons_tmp_unsorted_s(i).coords(:,2),neurons_tmp_unsorted_s(i).coords(:,1),'color',colormap_s(i,:));
        end
        axis equal
        while check && go==0
            drawnow
        end
    end
    save([outputpath '\parameters3.mat'],'max_gap_closing_distance','max_interframe_distance','non_linking_cost','neurons_tmp_unsorted_s','-v7.3');
else
    load([outputpath '\parameters3.mat']);
end
try
    structfun(@delete,guiapp);
    delete(text_non_linking_cost);
    delete(text_max_gap_closing_distance);
    delete(text_max_interframe_distance);
    delete(text_length_filter);
end

%% Step 2.1: creating distance matrix of subvideo tracks

check=0;
go=0;
try
    clear global short_neurons_tmp_unsorted_s
end
global neurons_bv
global short_neurons_tmp_unsorted_s
for i=1:numel(neurons_tmp_unsorted_s)
    x_vector=neurons_tmp_unsorted_s(i).coords(:,1);
    non_nan_times=find(~isnan(x_vector));
    short_neurons_tmp_unsorted_s(i).coords(:,:)=neurons_tmp_unsorted_s(i).coords(non_nan_times,:);
end
clear distmatrix_s
distmatrix_s=makeDistMatrix_s(neurons_tmp_unsorted_s,maxtime_s);

%% Distance matrix parameter selection in subvideo
try
    delete(guiapp.non_linking_cost)
    delete(guiapp.max_gap_closing_distance)
    delete(guiapp.max_interframe_distance)
    delete(text_non_linking_cost)
    delete(text_max_gap_closing_distance)
    delete(text_max_interframe_distance)
    delete(text_length_filter)
end

max_gap_closing_mutual_distance=4;
global go
global check
global display_matrix
go=0;
check=0;
framenum(framenum>DURATA_TOTALE_s)=DURATA_TOTALE_s;
if ~isfile([outputpath '\parameters4.mat'])
    while check==0 && go==0 % qui inserire il salvataggio dei dati
        check=1;
        try
            figure(f)
        catch
            f=figure('Units','normalized','Position',[0.125 0.125 0.75 0.75]);
        end
        % Step 2.2: gap closing based on mutual distance
        % parameters
        mutual_position_based_links=[]; reconstructed_links=[];
        mutual_position_based_links=get_mpbl(distmatrix_s,short_neurons_tmp_unsorted_s,maxtime_s,xy_pixel_size,z_pixel_size,max_gap_closing_mutual_distance);
        reconstructed_links=reconstruct_links(mutual_position_based_links,max_gap_closing_mutual_distance);
        neurons_bv=create_neurons_bv(short_neurons_tmp_unsorted_s,reconstructed_links);

        try
            shown_fig.CData=I;
        catch
            global shown_fig
            global ax
            shown_fig=imagesc(I);
            ax=gca;
            axis equal
            axis off
        end
        % slider range definition
        maxsize=25; minsize=2; sizerange=maxsize-minsize;
        minsigma=3; maxsigma=150; sigmarange=maxsigma-minsigma;
        mintime=1; timerange=maxtime_s-mintime;
        minslice=1; maxslice=number_of_planes; slicerange=maxslice-minslice; slicerange(slicerange==0)=1;
        disp('Tracking parameters selection...')
        guiapp.sld_timepoint= uicontrol('parent',f,'Style','slider','Value',framenum,'min',mintime,'max',maxtime_s,'SliderStep',[1 5]./timerange,'Units','normalized','Position',[0.05 0.1525 0.9 0.025]);
        guiapp.sld_slice= uicontrol('parent',f,'Style','slider','Value',slice,'min',minslice,'max',maxslice,'SliderStep',[1 5]./slicerange,'Units','normalized','Position',[0.05 0.1275 0.9 0.025]);
        text_timepoint = uicontrol('parent',f,'Style','text','string',['framenum = ' num2str(framenum)],'Position',[650 35 150 25]);
        text_slice = uicontrol('parent',f,'Style','text','string',['slice = ' num2str(slice)],'Position',[650 10 150 25]);
        text_spot_num = uicontrol('parent',f,'Style','text','string',['spot number = ' num2str(0)],'Position',[650 60 150 25]);
        set(guiapp.sld_timepoint,'Callback','global framenum; framenum=round(guiapp.sld_timepoint.Value); sliderMoving_tracks_s_3;');
        set(guiapp.sld_slice,'Callback','global slice; slice=guiapp.sld_slice.Value; sliderMoving_tracks_s_3;');
        guiapp.max_gap_closing_mutual_distance = uicontrol('parent',f,'Style','edit','String',num2str(max_gap_closing_mutual_distance),'Position',[150 20 150 25]);
        text_max_gap_closing_mutual_distance = uicontrol('parent',f,'Style','text','string',['max gap closing mutual distance = ' num2str(max_gap_closing_mutual_distance)],'Position',[150 55 150 25]);
        set(guiapp.max_gap_closing_mutual_distance,'Callback','global max_gap_closing_mutual_distance; max_gap_closing_mutual_distance=str2num(guiapp.max_gap_closing_mutual_distance.String); set(text_max_gap_closing_mutual_distance,"string",["max gap closing mutual distance = " (guiapp.max_gap_closing_mutual_distance.String)]);');
        guiapp.bttn_go= uicontrol('parent',f,'Style','pushbutton','String','GO!','Units','normalized','Position',[0.95 0.0475 0.04 0.05],'Callback','global go; go=1;');
        guiapp.bttn_check= uicontrol('parent',f,'Style','pushbutton','String','Check','Value',0,'Units','normalized','Position',[0.85 0.0475 0.04 0.05],'Callback','global check; check=0;');
        colormap('gray')
        colormap_s=distinguishable_colors(numel(neurons_bv));
        figure(f)
        try
            delete(plotted_tracks);
        end
        plotted_tracks=[];
        hold on
        for i=1:numel(neurons_bv)
            plotted_tracks(i)=plot(neurons_bv(i).coords(:,2),neurons_bv(i).coords(:,1),'color',colormap_s(i,:));
            hold on
        end
        display_matrix=[];
        display_matrix=NaN.*ones(numel(neurons_bv),maxtime_s,3);
        for i=1:numel(neurons_bv)
            display_matrix(i,neurons_bv(i).coords(:,4),1:3)=neurons_bv(i).coords(:,1:3);
        end
        sliderMoving_tracks_s_3
    end
    save([outputpath '\parameters4.mat'],'thresholds','sizes','sigmas','zsizes','time_processing','max_gap_closing_mutual_distance','-v7.3');
else
    load([outputpath '\parameters4.mat']);
end
close all

