%% Initialization

for initialization=1
    addpath('C:\Users\enric\OneDrive\Documents\SeeElegans\Functions')

    % initial parameters
    global xy_pixel_size; global z_pixel_size
    % parameters to be chosen
    global thresholds; global sizes; global zsizes; global sigmas; global time_processing
    % visualization
    global framenum; global I; global slice; global ax
    global stackpath; global wholestack; global DURATA_TOTALE
    global NUMERO_PIANI
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
    stackpath='C:\Users\enric\OneDrive\Documents\SeeElegans\Dataset1_Test\';
    outputpath=['C:\Users\enric\OneDrive\Documents\SeeElegans\Dataset1_Test\output\'];

    if ~exist(outputpath, 'dir')
        mkdir(outputpath)
    end


    load([stackpath 'chosen_ws.mat']);
    if size(wholestack,1)>size(wholestack,2)
        wholestack=double(permute(wholestack,[2 1 3 4]));
    end

    NUMERO_PIANI=size(wholestack,3);
    DURATA_TOTALE=size(wholestack,4);
    maxtime=DURATA_TOTALE;

    audio=0;
end

%% Step 0.0: visualization and segmentation parameter selection

for collapsing=1
    global go;

    for collapse=1
        f=figure('Units','normalized','Position',[0.125 0.125 0.75 0.75]);
        set(f, 'MenuBar', 'none');
        ax=axes(f,'Position',[0.05 0.2 0.9 0.8]);
        framenum(framenum>DURATA_TOTALE)=DURATA_TOTALE;
        I=wholestack(:,:,slice,framenum)-mean(mean(wholestack(1:38,1:38,slice,framenum)));
        imagesc(I), axis equal, axis off

        % slider range definition
        max_threshold=max(wholestack(:));
        maxsize=25; minsize=2; sizerange=maxsize-minsize;
        minsigma=1; maxsigma=150; sigmarange=maxsigma-minsigma;
        mintime=1; maxtime=DURATA_TOTALE; timerange=maxtime-mintime;
        minslice=1; maxslice=NUMERO_PIANI; slicerange=maxslice-minslice; slicerange(slicerange==0)=1;

        disp('Starting...')

        guiapp.sld_threshold = uicontrol('parent',f,'Style','slider','Value',max_threshold./2,'min',2,'max',max_threshold,'Position',[50 10 150 25]);
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
    end


    if ~isfile([outputpath '\parameters.mat'])
        while go==0
            drawnow
        end
        save([outputpath '\parameters.mat'],'thresholds','sizes','sigmas','zsizes','time_processing','-v7.3');
    else
        load([outputpath '\parameters.mat'])
    end

    structfun(@delete,guiapp);
end

%% Step 0.1: visualization and subvideo cropping

for collapsing=1
    global parameters
    global texts
    global go;
    global check;

    parameters.first_frame=1;
    parameters.last_frame=maxtime;
    guiapp.first_frame = uicontrol('parent',f,'Style','pushbutton','Value',1,'Position',[1050 60 125 25],'String','Start Frame');
    guiapp.last_frame = uicontrol('parent',f,'Style','pushbutton','Value',maxtime,'Position',[1050 30 125 25],'String','End Frame');
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
    %         if ~isfile([outputpath '\parameters2.mat'])
    %             while go==0
    %                 drawnow
    %             end
    %             disp(['Evaluating from ' num2str(parameters.first_frame) ' to ' num2str(parameters.last_frame) '...']);
    %             wholestack_s=wholestack(:,:,:,parameters.first_frame:parameters.last_frame);
    %             save([outputpath '\parameters2.mat'],'parameters','-v7.3');
    %         else
    %             load([outputpath '\parameters2.mat']);
    %
    %         end
    wholestack_s=wholestack(:,:,:,parameters.first_frame:parameters.last_frame);
    structfun(@delete,guiapp);
    structfun(@delete,texts);
end

%% Step 1: finding spots (LoG filtering and thresholding) in subvideo

for collapsing=1
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
        % finding spots ----------- MEGLIO FARLA QUI LA SUBPIXEL LOCALIZATION
        disp('Step 1.1: finding spots')
        spots_s{local_maxtime}=[];

        global p, global h
        D = parallel.pool.DataQueue;
        h = waitbar(0, 'Finding spots ...');
        p=1;
        afterEach(D, @nUpdateWaitbar_s);
        for framenum_id=1:local_maxtime
            spots_s{framenum_id}=bd3d_t_local_s(wholestack_s(:,:,:,framenum_id),sizes_local,thresholds_local,sigmas_local,zsizes_local);
            send(D, framenum_id);
        end
        elapsed_time_step_1=toc;
        A_time_s=spots_s;
        disp(['computation time: ' num2str(elapsed_time_step_1)]);
        close(h)
        clear global p
        clear global h
    else
        if ~exist('A_time_s','var')
            load([outputpath '\A_time_s.mat']);
        end
    end
end

clear neurons_reconstructed_max_2
for i=1:numel(A_time_s{1,1})
    neurons_reconstructed_max_2(i).coords=A_time_s{1,1}(i).center;
end
save([outputpath '\neurons_reconstructed_max_2.mat'],'neurons_reconstructed_max_2','thresholds','sizes','zsizes','sigmas','elapsed_time_step_1','-v7.3');
