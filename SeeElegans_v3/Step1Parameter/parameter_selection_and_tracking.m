function parameter_selection_and_tracking(data, outputpath, voxel_size)

% tasto deseleziona tutto
% tasto per la scala dei colori nel video
% tasto luminosità
% deve lasciare un simbolo nel punto in cui tracki
% se i primi vicini sono troppo lontani non deve considerarli
% se passi col mouse su una traccia sarebbe carino vederla nel tempo

% MANCA ZSIZES LOCAL, LO METTO O NO?
% NON HO MESSO I PARAMETRI DELLO STEP 1, PRIMO PANNELLO DIPENDENTI DALLA
% VOXEL SIZE. 

% Parameter initialization ------------------------------------------------
addpath('./Functions/');
stack_initialization();

% GUI
zmax=size(data,3);
tmax=size(data,4);
tmax_local=[];
z=1;
t=1;
menu_control=1;
if nargin<3 || numel(voxel_size)<3
    voxel_size=[0.267 0.267 2];
end
zratio=voxel_size(3)/voxel_size(1);
factor = 10;

% Parameters 1
thresholds=200;
zsizes=3;
sigmas=0.75;
sizes=2 * ceil(3 * sigmas) + 1;

% Parameters 2
data_s=[];
parameters.first_frame=1;
parameters.last_frame=round(tmax/10);

% Parameters 3
spots_s{1}=[];
non_linking_cost=1;
max_gap_closing_distance=1;
max_interframe_distance=2;
length_filter=1;
coords_t=[];
coords_colors=[];

% Parameters 4
coords_t2=[];
coords_colors_2=[];
neurons_bv=[];
max_gap_closing_mutual_distance=2;

% Create a figure window
fig = figure('visible','off');

FontSize_factor=2;
current_size=get(fig, 'position');
guiapp.control_panel=uipanel('Parent',fig,'Units','normalized','Position',[0 0.2 1 0.15],'BorderType','none');

% Create axes to display image
guiapp.img_ax = axes('Parent', fig, 'Position', [0.1 0.3 0.8 0.6]);
guiapp.img_slice = imagesc(data(:,:,z,t));
set(guiapp.img_ax,'Units','pixels')
currently_displayed_pxsz=diff(guiapp.img_ax.XLim)/guiapp.img_ax.Position(3);
set(guiapp.img_ax,'Units','normalized')

hold(guiapp.img_ax,'on');
scale_shift=5;
displayed_scale=10;
scl_x=[0 displayed_scale/voxel_size(1)]+scale_shift;
scl_y=[0 0]+scale_shift;
text(mean(scl_x),mean(scl_y)+8,[num2str(displayed_scale) ' \mum'],'Color','w','Interpreter','tex','HorizontalAlignment','center','FontSize',9);
guiapp.scalebar = plot(scl_x,scl_y,'Color','w','LineWidth',2); 
hold(guiapp.img_ax,'on');
guiapp.scatter_plot = scatter(NaN,NaN,NaN);
hold(guiapp.img_ax,'on');
guiapp.track_plot = scatter(NaN,NaN,NaN);
update_spots=@update_spots_step_1_and_2;
update_image=@update_image_t;
axis equal
axis off
colormap(gray);

% Create sliders for changing the Z and t dimensions
guiapp.zSlider = uicontrol('Parent', fig, 'units', 'normalized', 'Style', 'slider', 'Min', 1, 'Max', zmax, 'Value', z, ...
    'SliderStep', [1/(zmax-1) 5/(zmax-1)], 'Position', [0.2 0.01 0.6 0.04], 'Callback', @zSlider_callback);
guiapp.tSlider = uicontrol('Parent', fig, 'units', 'normalized', 'Style', 'slider', 'Min', 1, 'Max', tmax, 'Value', t, ...
    'SliderStep', [1/(tmax-1) 5/(tmax-1)], 'Position', [0.2 0.05 0.6 0.04], 'Callback', @tSlider_callback);

% Create text boxes to display the current Z, t, and number of spot
guiapp.zValue = uicontrol('Parent', fig, 'Style', 'text', 'String', ['z: ' num2str(get(guiapp.zSlider, 'Value'))], ...
    'units', 'normalized', 'Position', [0.8 0.01 0.05 0.04],'HorizontalAlignment','left');
guiapp.tValue = uicontrol('Parent', fig, 'Style', 'text', 'String', ['t: ' num2str(get(guiapp.tSlider, 'Value'))], ...
    'units', 'normalized', 'Position', [0.8 0.05 0.08 0.04],'HorizontalAlignment','left');
% Advance
guiapp.advance_button = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', char(8594), ... % right arrow
    'Units', 'normalized', 'Position', [0.94 0.0 0.06 0.05], 'visible', 'off', 'Callback', @next_menu);
guiapp.goback_button = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', char(8592), ... % right arrow
    'Units', 'normalized', 'Position', [0.88 0.0 0.06 0.05], 'visible', 'off', 'Callback', @previous_menu);
% Save
guiapp.save_parameters = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'save', ...
    'Units', 'normalized', 'Position', [0.88 0.05 0.12 0.05]);

set(fig, 'visible','on');
set(fig, 'WindowScrollWheelFcn', @mouse_scroll_callback_t);
set(fig, 'WindowKeyPressFcn', @keypress_callback);
set(fig, 'WindowKeyReleaseFcn', @keyrelease_callback);
set(fig, 'SizeChangedFcn', @resize_text_function);
addlistener(guiapp.img_ax,'XLim','PostSet',@resize_spots_function);

% Set Menu COntrols--------------------------------------------------------
setmenu();
resize_text_function();

% Functions ---------------------------------------------------------------
    function stack_initialization
        data=single(data);
        if size(data,1)>size(data,2)
            data=single(permute(data,[2 1 3 4]));
        end
    end

    % Mouse behaviour
    
    function keypress_callback(~, event)
        if strcmp(event.Key, 'control')
            set(fig, 'WindowScrollWheelFcn', @mouse_scroll_callback_z);
        end
    end
    
    function keyrelease_callback(~, event)
        if strcmp(event.Key, 'control')
            set(fig, 'WindowScrollWheelFcn', @mouse_scroll_callback_t);
        end
    end
    
    function mouse_scroll_callback_z(~, event)
        z = z + event.VerticalScrollCount;
        z = max(min(z, zmax), 1);
        guiapp.zSlider.Value = z;
        guiapp.zValue.String = sprintf('z: %d', z);
        update_image();
    end
    
    function mouse_scroll_callback_t(~, event)
        t = t + round(10.*event.VerticalScrollCount);
        t = max(min(t, guiapp.tSlider.Max), 1);
        guiapp.tSlider.Value = t;
        guiapp.tValue.String = sprintf('t: %d', t);
        update_image();
    end

    
    % Callback functions for the buttons and sliders
    
    function zSlider_callback(source, ~)
        z = round(get(source, 'Value'));
        guiapp.zValue.String = sprintf('z: %d', z);
        update_image();
    end
    
    function tSlider_callback(source, ~)
        t = round(get(source, 'Value'));
        guiapp.tValue.String = sprintf('t: %d', t);
        update_image();
    end

    % Callback functions for resizing
    
    function resize_text_function(~,~)
        new_size=get(fig,'Position');
        [~, changed_dim]=min(new_size(3:4));
        FontSize_factor=10*new_size(2+changed_dim)./current_size(2+changed_dim);
        textObjs = findobj(fig, '-property', 'FontSize');
        for i=1:numel(textObjs)
            set(textObjs(i),'FontSize',FontSize_factor);
        end
        resize_spots_function();
    end

    function resize_spots_function(~,~)
        set(guiapp.img_ax,'Units','pixels')
        new_displayed_pxsz=guiapp.img_ax.Position(3)/diff(guiapp.img_ax.XLim);
        set(guiapp.img_ax,'Units','normalized');
        factor=factor*new_displayed_pxsz/(currently_displayed_pxsz);
        try
            update_spots();
        end
        currently_displayed_pxsz=new_displayed_pxsz;
    end

    % Functions for data display

    function update_image_t
        set(guiapp.img_slice,'CData',(data(:,:,z,t)));
        drawnow
        update_spots();
    end

    function update_image_s
        set(guiapp.img_slice,'CData',(data_s(:,:,z,t)));
        drawnow
        update_spots();
    end

    function update_spots_step_1_and_2(~,~)
        delete(guiapp.scatter_plot)
        if guiapp.step(1).show_spots_check.Value==1
            spots=bd3d(data(:,:,:,t),thresholds,sizes,zsizes,sigmas,voxel_size);
            guiapp.scatter_plot=plotspots(guiapp.img_ax, guiapp.scatter_plot, spots, sizes, z, zratio, factor);
            guiapp.step(1).set_spotnum_text.String=sprintf('%d',size(spots,1));
        end
    end

    function update_spots_step_3(~,~)
        delete(guiapp.scatter_plot)
        current_coords=squeeze(coords_t(t,:,:))';
        guiapp.scatter_plot=plotspots(guiapp.img_ax, guiapp.scatter_plot, current_coords, sizes, z, zratio, factor, coords_colors);
        zmarker=abs(squeeze(coords_t(t,3,:))-z);
        set(guiapp.track_plot,'LineWidth',0.8);
        set(guiapp.track_plot(zmarker<2),'LineWidth',1.2);
        set(guiapp.track_plot(zmarker<1),'LineWidth',2.2);

        guiapp.step(3).set_spotnum_text.String=sprintf('%d',size(spots_s{t},1));
    end

    function update_spots_step_4(~,~)
        delete(guiapp.scatter_plot)
        current_coords=squeeze(coords_t2(t,:,:))';
        guiapp.scatter_plot=plotspots(guiapp.img_ax, guiapp.scatter_plot, current_coords, sizes, z, zratio, factor, coords_colors_2);
        zmarker=abs(squeeze(coords_t2(t,3,:))-z);
        set(guiapp.track_plot,'LineWidth',0.8);
        set(guiapp.track_plot(zmarker<2),'LineWidth',1.2);
        set(guiapp.track_plot(zmarker<1),'LineWidth',2.2);

        guiapp.step(4).set_spotnum_text.String=sprintf('%d',size(spots_s{t},1));
    end

    
    % Functions to control behaviour of menu

    function overwrite_answer=overwrite_request(overwriting_path)
        overwriting_path=strrep(overwriting_path,'\','/');
        overwriting_path=strrep(overwriting_path,'//','/');
        overwrite_answer=questdlg(['Do you want to overwrite existing parameter settings in the following path?' newline overwriting_path ]);
        switch overwrite_answer
            case 'Yes'
                overwrite_answer=1;
        end
    end

    function save_parameters_step_1(~,~)
        if ~exist(outputpath,'dir')
            mkdir(outputpath)
        end
        
        if ~exist([outputpath '/parameters.mat'],"file")
            save([outputpath '/parameters.mat'],'thresholds','sizes','sigmas','zsizes','voxel_size','-v7.3');
            set(guiapp.advance_button,'visible','on');
        else
            overwrite_answer=overwrite_request([outputpath '/parameters.mat']);
            if overwrite_answer==1
                save([outputpath '/parameters.mat'],'thresholds','sizes','sigmas','zsizes','voxel_size','-v7.3');
            end
        end
    end

    function save_parameters_step_2(~,~)
        if ~exist([outputpath '/parameters2.mat'],"file")
            save([outputpath '/parameters2.mat'],'parameters','-v7.3');
            set(guiapp.advance_button,'visible','on')
        else
            overwrite_answer=overwrite_request([outputpath '/parameters2.mat']);
            if overwrite_answer==1
                save([outputpath '/parameters2.mat'],'parameters','-v7.3');
                data_s=[];
                data_s=data(:,:,:,parameters.first_frame:parameters.last_frame);
            end
        end
    end

    function save_parameters_step_3(~,~)
        if ~exist([outputpath '/parameters3.mat'],"file")
            save([outputpath '/parameters3.mat'],'non_linking_cost','max_gap_closing_distance','max_interframe_distance','length_filter','-v7.3');
            set(guiapp.advance_button,'visible','on')
        else
            overwrite_answer=overwrite_request([outputpath '/parameters3.mat']);
            if overwrite_answer==1
                save([outputpath '/parameters3.mat'],'non_linking_cost','max_gap_closing_distance','max_interframe_distance','length_filter','-v7.3');
            end
        end
    end

    function save_parameters_step_4(~,~)
        if ~exist([outputpath '/parameters4.mat'],"file")
            save([outputpath '/parameters4.mat'],'max_gap_closing_mutual_distance',"neurons_bv",'-v7.3');
            set(guiapp.advance_button,'visible','on')
        else
            overwrite_answer=overwrite_request([outputpath '/parameters4.mat']);
            if overwrite_answer==1
                save([outputpath '/parameters4.mat'],'max_gap_closing_mutual_distance',"neurons_bv",'-v7.3');
            end
        end
    end

    function next_menu(~,~)
        menu_control=menu_control+1;
        menu_control(menu_control > 5 )=5;
        setmenu
    end

    function previous_menu(~,~)
        menu_control=menu_control-1;
        setmenu
    end

    function setmenu(~,~)
        switch menu_control
            case 1
                delete(guiapp.control_panel.Children)
                guiapp.tSlider.Value=round(size(data,4)/2);
                guiapp.tSlider.Max=size(data,4);
                guiapp.tSlider.SliderStep=[1/(tmax-1) 5/(tmax-1)];
                update_image=@update_image_t;
                step1_GUI();
                guiapp.save_parameters.Callback=@save_parameters_step_1;
                set(guiapp.goback_button,'visible','off')
                if exist([outputpath '/parameters.mat'],"file")
                    set(guiapp.advance_button,'visible','on')
                end
                update_spots();
                resize_text_function();

            case 2
                delete(guiapp.control_panel.Children)
                guiapp.tSlider.Value=round(size(data,4)/2);
                guiapp.tSlider.Max=size(data,4);
                guiapp.tSlider.SliderStep=[1/(tmax-1) 5/(tmax-1)];
                update_image=@update_image_t;
                step2_GUI();
                guiapp.save_parameters.Callback=@save_parameters_step_2;
                set(guiapp.goback_button,'visible','on')
                if exist([outputpath '/parameters2.mat'],"file")
                    set(guiapp.advance_button,'visible','on')
                end
                update_spots();
                resize_text_function();

            case 3
                set(guiapp.advance_button,'visible','off')
                if exist([outputpath '/parameters3.mat'],"file")
                    set(guiapp.advance_button,'visible','on')
                    load([outputpath '/parameters3.mat'],'non_linking_cost','max_gap_closing_distance','max_interframe_distance','length_filter');
                end
                update_image=@update_image_s;
                initialize_step_3();
                delete(guiapp.control_panel.Children)
                guiapp.tSlider.Value=round(size(data_s,4)/2);
                t=guiapp.tSlider.Value;
                guiapp.tSlider.Max=size(data_s,4);
                guiapp.tSlider.SliderStep=[1/(tmax_local-1) 5/(tmax_local-1)];
                step3_GUI();
                guiapp.save_parameters.Callback=@save_parameters_step_3;
                set(guiapp.goback_button,'visible','on')
                generate_tracks();
                update_spots=@update_spots_step_3;
                update_spots();
                update_image();
                resize_text_function();

            case 4
                delete(guiapp.control_panel.Children)
                guiapp.tSlider.Value=round(size(data_s,4)/2);
                guiapp.tSlider.Max=size(data_s,4);
                guiapp.tSlider.SliderStep=[1/(tmax_local-1) 5/(tmax_local-1)];
                update_image=@update_image_s;
                step4_GUI();
                if isfile([outputpath '/parameters4.mat'])
                    load([outputpath '/parameters4.mat'],'max_gap_closing_mutual_distance');
                    guiapp.step(4).max_gap_closing_mutual_distance.String=num2str(max_gap_closing_mutual_distance);
                    set(guiapp.advance_button,'visible','on')
                end
                guiapp.save_parameters.Callback=@save_parameters_step_4;
                set(guiapp.goback_button,'visible','on')
                regenerate_tracks();
                update_spots=@update_spots_step_4;
                update_spots();
                resize_text_function();

            case 5
                disp('Tracking Started')
                delete(guiapp.control_panel.Children)
                set(guiapp.goback_button,'visible','off')
                set(guiapp.advance_button,'visible','off')
                msgbox('Tracking Started')
                drawnow
                tracking_step(data, outputpath);
        end
    end

    function step1_GUI
        % thresholds
        guiapp.step(1).set_threshold = uicontrol('Parent', guiapp.control_panel, 'Style', 'edit', 'String', num2str(thresholds), ...
            'Units', 'normalized', 'Position', [0.08 0.2 0.12 0.4], 'Callback', @button_thresholds);
        guiapp.step(1).set_threshold_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Threshold', ...
            'Units', 'normalized', 'Position', [0.08 0.6 0.12 0.4], 'Callback', @button_thresholds);
        % sizes
        guiapp.step(1).set_size = uicontrol('Parent', guiapp.control_panel, 'Style', 'edit', 'String', num2str(sizes), ...
            'Units', 'normalized', 'Position', [0.2 0.2 0.12 0.4], 'Callback', @button_sizes);
        guiapp.step(1).set_size_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Size', ...
            'Units', 'normalized', 'Position', [0.2 0.6 0.12 0.4], 'Callback', @button_thresholds);
        % sigmas
        guiapp.step(1).set_sigma = uicontrol('Parent', guiapp.control_panel, 'Style', 'edit', 'String', num2str(sigmas), ...
            'Units', 'normalized', 'Position', [0.32 0.2 0.12 0.4], 'Callback', @button_sigmas);
        guiapp.step(1).set_sigma_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Sigma', ...
            'Units', 'normalized', 'Position', [0.32 0.6 0.12 0.4]);
        % number of spots
        guiapp.spotnum_title = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Spot Num', ...
            'Units', 'normalized', 'Position', [0.44 0.6 0.12 0.4], 'Callback', @button_thresholds);
        guiapp.step(1).set_spotnum_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Spot Num', ...
            'Units', 'normalized', 'Position', [0.44 0.2 0.12 0.4]);
        guiapp.step(1).show_spots_check = uicontrol('Parent', guiapp.control_panel, 'Style', 'checkbox', 'String', 'spots', ...
            'Units', 'normalized', 'Position', [0 0 0.09 0.25], 'Value', 1, 'Callback', @check_spots);
        set(guiapp.advance_button,'visible','off')
        
        function button_thresholds(source, ~)
            thresholds = str2double(source.String);
            update_spots();
        end

        function button_sizes(source, ~)
            sizes = str2double(source.String);
            update_spots();
        end

        function button_sigmas(source, ~)
            sigmas = str2double(source.String);
            suggested_size = 2 * ceil(3 * sigmas) + 1;
            guiapp.step(1).set_size.String = num2str(suggested_size);
            sizes = suggested_size;
            update_spots();
        end

        function check_spots(~,~)
            update_spots();
        end
        
        update_spots=@update_spots_step_1_and_2;
    end

    function step2_GUI
        guiapp.step(2).first_frame = uicontrol('parent', guiapp.control_panel,'Style','pushbutton','Units', 'normalized','Position',[0.08 0.2 0.12 0.4],'String','Set Start','Callback',@setStartFrame);
        guiapp.step(2).last_frame = uicontrol('parent', guiapp.control_panel,'Style','pushbutton','Units', 'normalized','Position',[0.2 0.2 0.12 0.4],'String','Set Stop','Callback',@setEndFrame);
        guiapp.step(2).first_frame_texts = uicontrol('parent', guiapp.control_panel,'Style','text','string',['Start = ' num2str(parameters.first_frame)],'Units', 'normalized','Position',[0.08 0.6 0.12 0.4]);
        guiapp.step(2).last_frame_texts = uicontrol('parent', guiapp.control_panel,'Style','text','string',['Stop = ' num2str(parameters.last_frame)],'Units', 'normalized','Position',[0.2 0.6 0.12 0.4]);
        
        % number of spots
        guiapp.step(2).spotnum_title = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Spot Num', ...
            'Units', 'normalized', 'Position', [0.44 0.6 0.12 0.4], 'Callback', @button_thresholds);
        guiapp.step(1).set_spotnum_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', '', ...
            'Units', 'normalized', 'Position', [0.44 0.2 0.12 0.4], 'Callback', @button_thresholds);
        set(guiapp.advance_button,'visible','off')
        guiapp.step(1).show_spots_check = uicontrol('Parent', guiapp.control_panel, 'Style', 'checkbox', 'String', 'spots', ...
            'Units', 'normalized', 'Position', [0 0 0.09 0.25], 'Value', 1, 'Callback', @check_spots);

        function setStartFrame(~,~)
            parameters.first_frame=t;
            guiapp.step(2).first_frame_texts.String=['Start = ' num2str(parameters.first_frame)];
        end

        function setEndFrame(~,~)
            parameters.last_frame=t;
            guiapp.step(2).last_frame_texts.String=['Stop = ' num2str(parameters.last_frame)];
        end
        
        function check_spots(~,~)
            update_spots();
        end

        update_spots=@update_spots_step_1_and_2;
    end

    function step3_GUI
        guiapp.step(3).non_linking_cost = uicontrol('parent', guiapp.control_panel,'Style','edit','String',num2str(non_linking_cost),'units','normalized','Position',[0.08 0.2 0.12 0.4],'Callback',@non_linking_cost_btn);
        guiapp.step(3).max_gap_closing_distance = uicontrol('parent', guiapp.control_panel,'Style','edit','String',num2str(max_gap_closing_distance),'units','normalized','Position',[0.2 0.2 0.12 0.4],'Callback',@max_gap_closing_distance_btn);
        guiapp.step(3).max_interframe_distance = uicontrol('parent', guiapp.control_panel,'Style','edit','String',num2str(max_interframe_distance),'units','normalized','Position',[0.32 0.2 0.12 0.4],'Callback',@max_interframe_distance_btn);
        guiapp.step(3).length_filter = uicontrol('parent', guiapp.control_panel,'Style','edit','String',num2str(length_filter),'units','normalized','Position',[0.44 0.2 0.12 0.4],'Callback',@length_filter_btn);        
        guiapp.step(3).check_button = uicontrol('parent', guiapp.control_panel,'Style','pushbutton','String','Check','units','normalized','Position',[0.84 0.2 0.12 0.4],'Callback',@generate_tracks);

        guiapp.step(3).text_non_linking_cost = uicontrol('parent', guiapp.control_panel,'Style','text','string',['Non linking cost = ' num2str(non_linking_cost)],'units','normalized','Position',[0.08 0.6 0.12 0.4]);
        guiapp.step(3).text_max_gap_closing_distance = uicontrol('parent', guiapp.control_panel,'Style','text','string',['Gap closing dist = ' num2str(max_gap_closing_distance)],'units','normalized','Position',[0.2 0.6 0.12 0.4]);
        guiapp.step(3).text_max_interframe_distance = uicontrol('parent', guiapp.control_panel,'Style','text','string',['Inter-frame dist = ' num2str(max_interframe_distance)],'units','normalized','Position',[0.32 0.6 0.12 0.4]);
        guiapp.step(3).text_length_filter = uicontrol('parent', guiapp.control_panel,'Style','text','string',['Length filter = ' num2str(length_filter)],'units','normalized','Position',[0.44 0.6 0.12 0.4]);

        guiapp.step(3).spotnum_title = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Spot Num', 'Units', 'normalized', 'Position', [0.56 0.6 0.12 0.4]);
        guiapp.step(3).set_spotnum_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', '', 'Units', 'normalized', 'Position', [0.56 0.2 0.12 0.4]);
        guiapp.step(3).tracknum_title = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Track Num', 'Units', 'normalized', 'Position', [0.68 0.6 0.12 0.4]);
        guiapp.step(3).set_tracknum_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', '', 'Units', 'normalized', 'Position', [0.68 0.2 0.12 0.4]);

        function length_filter_btn(~,~)
            length_filter=str2double(guiapp.step(3).length_filter.String);
            guiapp.step(3).text_length_filter.String=['Length filter = ' num2str(length_filter)];
        end
        function max_interframe_distance_btn(~,~)
            max_interframe_distance=str2double(guiapp.step(3).max_interframe_distance.String);
            guiapp.step(3).text_max_interframe_distance.String=['Inter-frame dist = ' num2str(max_interframe_distance)];
        end
        function max_gap_closing_distance_btn(~,~)
            max_gap_closing_distance=str2double(guiapp.step(3).max_gap_closing_distance.String);
            guiapp.step(3).text_max_gap_closing_distance.String=['Gap closing dist = ' num2str(max_gap_closing_distance)];
        end
        function non_linking_cost_btn(~,~)
            non_linking_cost=str2double(guiapp.step(3).non_linking_cost.String);
            guiapp.step(3).text_non_linking_cost.String=['Non linking cost = ' num2str(non_linking_cost)];
        end
    end

    function step4_GUI
        guiapp.step(4).max_gap_closing_mutual_distance = uicontrol('parent', guiapp.control_panel,'Style','edit','String',num2str(max_gap_closing_mutual_distance),'units','normalized','Position',[0.08 0.2 0.12 0.4],'Callback',@max_gap_closing_mutual_distance_btn);
        guiapp.step(4).check_button = uicontrol('parent', guiapp.control_panel,'Style','pushbutton','String','Check','units','normalized','Position',[0.84 0.2 0.12 0.4],'Callback',@regenerate_tracks);
        
        guiapp.step(4).texts_max_gap_closing_mutual_distance = uicontrol('parent', guiapp.control_panel,'Style','text','string',['gap closing dist = ' num2str(max_gap_closing_mutual_distance)],'Units', 'normalized','Position',[0.08 0.6 0.12 0.4]);
        
        guiapp.step(4).spotnum_title = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Spot Num', 'Units', 'normalized', 'Position', [0.44 0.6 0.12 0.4]);
        guiapp.step(4).set_spotnum_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', '', 'Units', 'normalized', 'Position', [0.44 0.2 0.12 0.4]);
        guiapp.step(4).tracknum_title = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', 'Track Num', 'Units', 'normalized', 'Position', [0.68 0.6 0.12 0.4]);
        guiapp.step(4).set_tracknum_text = uicontrol('Parent', guiapp.control_panel, 'Style', 'text', 'String', '', 'Units', 'normalized', 'Position', [0.68 0.2 0.12 0.4]);


        set(guiapp.advance_button,'visible','off')

        update_spots=@update_spots_step_4;

        function max_gap_closing_mutual_distance_btn(~,~)
            max_gap_closing_mutual_distance=str2double(get(guiapp.step(4).max_gap_closing_mutual_distance,'String'));
            guiapp.step(4).texts_max_gap_closing_mutual_distance.String=['gap closing dist = ' num2str(max_gap_closing_mutual_distance)];
        end
    end

    
    % Functions to process data

    function initialize_step_3(~,~)
        thresholds_local=thresholds;
        sigmas_local=sigmas;
        sizes_local=sizes;
        zsizes_local=zsizes;
        data_s=[];
        data_s=data(:,:,:,parameters.first_frame:parameters.last_frame);
        parameters_local=parameters;
        tmax_local=size(data_s,4);

        if ~isfile([outputpath '/spots_s.mat'])
            short_video_spots;
            save([outputpath '/spots_s.mat'],'spots_s','sizes_local','thresholds_local','sigmas_local','zsizes_local','parameters_local');
        else
            previous_settings=load([outputpath '/spots_s.mat']);
            if previous_settings.thresholds_local~=thresholds || ...
                    previous_settings.sigmas_local~=sigmas || ...
                    previous_settings.sizes_local~=sizes || ...
                    previous_settings.zsizes_local~=zsizes || ...
                    previous_settings.parameters_local.first_frame~=parameters.first_frame ||...
                    previous_settings.parameters_local.last_frame~=parameters.last_frame

                if overwrite_request([outputpath '/spots_s.mat'])==1
                    short_video_spots;
                    save([outputpath '/spots_s.mat'],'spots_s','sizes_local','thresholds_local','sigmas_local','zsizes_local','parameters_local');
                else
                    load_spots_from_file
                end
            else
                load_spots_from_file
            end
        end

        function load_spots_from_file
            disp('Loading from file');
            spots_s=previous_settings.spots_s;
            thresholds_local=previous_settings.thresholds_local;
            sigmas_local=previous_settings.sigmas_local;
            sizes_local=previous_settings.sizes_local;
            zsizes_local=previous_settings.zsizes_local;
            parameters.first_frame=previous_settings.parameters_local.first_frame;
            parameters.last_frame=previous_settings.parameters_local.last_frame;
            data_s=[];
            data_s=data(:,:,:,parameters.first_frame:parameters.last_frame);

        end

        function short_video_spots
            tic
            % finding spots -----------
            disp('Step 1.1: finding spots')
            spots_s{1}=[];
            h = waitbar(0, 'Finding spots ...');
            for framenum_id=1:tmax_local
                % spots_s{framenum_id}=bd3d_local(data_s(:,:,:,framenum_id),sizes_local,thresholds_local,sigmas_local,zsizes_local);
                spots_s{framenum_id}=bd3d(data_s(:,:,:,framenum_id),thresholds,sizes,zsizes,sigmas,voxel_size);
                waitbar(framenum_id/tmax_local,h);
            end
            close(h)
            elapsed_time_step_1=toc;
            disp(['computation time: ' num2str(elapsed_time_step_1)]);
        end

    end

    function generate_tracks(~,~)
        non_linking_cost = str2double(guiapp.step(3).non_linking_cost.String);
        max_gap_closing_distance = str2double(guiapp.step(3).max_gap_closing_distance.String);
        max_interframe_distance = str2double(guiapp.step(3).max_interframe_distance.String);
        length_filter = str2double(guiapp.step(3).length_filter.String);
        tmax_local=size(data_s,4);

        % LAP linker
        disp('Step 1.2: global minimum of interframe spatio-temporal spot proximity');
        h = waitbar(0, 'Interframe spot linking ...');
        spot_links_s{tmax_local}=[];
        tic
        for t_id=2:tmax_local
            spot_links_s{t_id}=LAP_linker_local_s(t_id,non_linking_cost,spots_s,voxel_size(1),voxel_size(3));
            waitbar(t_id/tmax_local,h);
        end
        close(h)
        toc

        % building segments lists
        disp('Step 1.3: building initial segments');
        [segmentlists_time_s,segmentlists_s]=build_links(spot_links_s);

        % LAP tracking: gap closing
        disp('Step 1.4: gap closing');
        segment_links_s=gap_closer_s(segmentlists_time_s,segmentlists_s,max_interframe_distance,max_gap_closing_distance,spots_s,voxel_size);

        % building track list
        tracklist_s=build_tracks(segment_links_s);

        % filter short linking lists
        disp('Step 1.5: trace length filter');
        filtered_exptracelist_s=merge_and_filter_short_links_s(tracklist_s,segmentlists_time_s,segmentlists_s,length_filter,tmax_local);

        % store tracks
        coords_colors=different_colors(numel(filtered_exptracelist_s));
        coords_t=NaN.*ones(tmax_local,4,numel(filtered_exptracelist_s));
        for i=1:numel(filtered_exptracelist_s)
            for t_id=1:tmax_local
                if ~isnan(filtered_exptracelist_s{i}(t_id)) && filtered_exptracelist_s{i}(t_id)~=0
                    id=filtered_exptracelist_s{i}(t_id);
                    sub_pixel_localization=spots_s{t_id}(id,1:3);
                    coords_t(t_id,:,i)=[sub_pixel_localization, t_id];
                end
            end
        end

        % update plots and text
        % colororder(guiapp.img_ax,different_colors(size(coords_t,3)));
        hold(guiapp.img_ax,'on');
        delete(guiapp.track_plot);
        guiapp.track_plot=plot(guiapp.img_ax,squeeze(coords_t(:,2,:)),squeeze(coords_t(:,1,:)),'Color','yellow');

        zmarker=abs(squeeze(coords_t(t,3,:))-z);

        set(guiapp.track_plot(zmarker<2),'LineWidth',2);
        set(guiapp.track_plot(zmarker<1),'LineWidth',4);

        set(guiapp.step(3).set_tracknum_text,'String',size(coords_t,3));
    end

    function regenerate_tracks(~,~)
        max_gap_closing_mutual_distance=str2double(get(guiapp.step(4).max_gap_closing_mutual_distance,'String'));
        guiapp.step(4).texts_max_gap_closing_mutual_distance.String=['gap closing dist = ' num2str(max_gap_closing_mutual_distance)];

        % Step 2.1: calculating mutual distance
        distmatrix_s=makeDistMatrix_s(coords_t);
        % Step 2.2: gap closing based on mutual distance parameters
        mutual_position_based_links=get_mpbl(distmatrix_s,coords_t,tmax_local,voxel_size(1),voxel_size(3),max_gap_closing_mutual_distance);
        reconstructed_links=reconstruct_links(mutual_position_based_links,max_gap_closing_mutual_distance);
        neurons_bv=create_neurons_bv(coords_t,reconstructed_links);
        coords_t2=nan(size(data_s,4),4,numel(neurons_bv));
        for i=1:numel(neurons_bv)
            valid_times=neurons_bv(i).coords(:,4);
            coords_t2(valid_times,:,i)=neurons_bv(i).coords;
        end

        % update plots and text
        coords_colors_2=different_colors(size(coords_t2,3));

        hold(guiapp.img_ax,'on');
        delete(guiapp.track_plot);
        guiapp.track_plot=plot(guiapp.img_ax,squeeze(coords_t2(:,2,:)),squeeze(coords_t2(:,1,:)),'Color','yellow');

        zmarker=abs(squeeze(coords_t2(t,3,:))-z);

        set(guiapp.track_plot(zmarker<2),'LineWidth',2);
        set(guiapp.track_plot(zmarker<1),'LineWidth',4);

        set(guiapp.step(4).set_tracknum_text,'String',size(coords_t2,3));
    end

end