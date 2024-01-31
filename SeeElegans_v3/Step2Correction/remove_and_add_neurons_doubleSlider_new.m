function remove_and_add_neurons(data, neurons_cleaning, outputpath, inversion)
% click on a neuron to delete it and then click on save to permanently 
% remove it from the sets of tracked spots
% click on button to add neuron and then save to check the resulting track

addpath('./Functions/');

% GUI Initialization
f = figure('Visible', 'off', 'Position', [360, 500, 450, 285]);
stack_initialization();

% variables initialization
guiapp=[];
timepoint = 1;
slice = 1;
selection=[];
to_remove=[];
to_add=[];
newly_tracked=[];
show_spots=1;
factor=100;
minca=min(data(:,:,:,1:20:end)); minca=min(minca(:));
maxca=max(data(:,:,:,1:20:end)); maxca=max(maxca(:));
custom_minca = minca;
custom_maxca = maxca;

initiate_GUI();

% plots initialization
hoover_over=plot(0,0,'.','Visible','off');
hoover_over_rt=plot(0,0,'.','Visible','off');
if nargin>3
    if inversion
        for n_id=1:numel(neurons_cleaning)
            neurons_cleaning(n_id).coords(:,[1 2])=neurons_cleaning(n_id).coords(:,[2 1]);
        end
    end
end
neurons_cleaning=pca_reconstruction(neurons_cleaning,1);
current_neuron_arangement= neurons_cleaning;
[X,Y,Z,C]=prepareCoordData;

% update the image display
image_ax=imagesc(data(:,:,slice,timepoint));
hold on
track_plot_ax=scatter(image_ax.Parent,X(:,timepoint),Y(:,timepoint),'or');
selected_track_plot_ax=scatter(image_ax.Parent,NaN,NaN);
hold on
colormap gray;
axis image;
axis off;

% trace plot
trace_ax=axes('Units','normalized','Position',[0.2 0.1 0.7 0.14],'YTickLabel',[],'XTickLabel',[]);
tm_xline=xline(trace_ax,timepoint);

axes(image_ax.Parent);
update_image();
set(f, 'Visible', 'on');

% Color Slider
hRangeSlider=caxis_slider();
color_scale_func();

% Mouse controls and Keypress
set(f, 'WindowScrollWheelFcn', @mouse_scroll_callback_t);
set(f, 'Windowkeypressfcn', @keypress_callback);
set(gcf, 'Keypressfcn', @keypress_callback,'KeyReleaseFcn',@keyrelease_callback);

% Resizing Commands
set(image_ax.Parent,'Units','pixels');
currently_displayed_pxsz=[diff(image_ax.Parent.XLim)/image_ax.Parent.Position(3) diff(image_ax.Parent.YLim)/image_ax.Parent.Position(4)];
set(image_ax.Parent,'Units','normalized');
FontSize_factor=[];
current_size=get(f, 'position');
set(f, 'SizeChangedFcn', @resize_text_function);
addlistener(image_ax.Parent,'XLim','PostSet',@resize_spots_function);
resize_text_function();

% functions ---------------------------------------------------------------
    % GUI initiation
    function initiate_GUI
        guiapp.sld_timepoint = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(data,4), 'Value', 1, 'Units', 'normalized', 'Position', [0.16, 0, 0.75, 0.04], 'Callback', @sliderMoving_timepoint);
        guiapp.sld_slice = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(data,3), 'SliderStep', [1/(size(data,3)-1) 4/(size(data,3)-1)], 'Value', 1, 'Units', 'normalized', 'Position', [0.16, 0.04, 0.75, 0.04], 'Callback', @sliderMoving_slice);
        guiapp.bttn_delneurons = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Del', 'FontSize', 6, 'Position', [0.005, 0.05, 0.1, 0.05], 'visible', 'off', 'Callback', @remove_neurons_func);
        guiapp.bttn_addneurons = uicontrol('Style', 'togglebutton', 'Units', 'normalized', 'String', 'Add', 'FontSize', 6, 'Position', [0.005, 0.1, 0.1, 0.05], 'visible', 'on', 'Callback', @update_controls);
        guiapp.bttn_addtracks = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Track', 'FontSize', 6, 'Position', [0.1125, 0.1, 0.1, 0.05], 'visible', 'off', 'Callback', @add_tracks_func);
        guiapp.bttn_savechanges = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Save', 'FontSize', 6, 'Position', [0.92, 0.01, 0.08, 0.05], 'visible', 'off', 'Callback', @save_neurons_func);
        guiapp.bttn_caxis_control = uicontrol('Style', 'checkbox', 'Value', 1, 'Units', 'normalized', 'String', 'scale', 'FontSize', 6, 'Position', [0.9, 0.92, 0.1, 0.05], 'visible', 'on', 'Callback', @color_scale_func);
        guiapp.bttn_spot = uicontrol('Style', 'checkbox', 'Value', 1, 'Units', 'normalized', 'String', 'Spots', 'FontSize', 6, 'Position', [0.01, 0.92, 0.1, 0.05], 'visible', 'on', 'Callback', @spot_func);
        guiapp.clear_selection = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Clear', 'FontSize', 6, 'Position', [0.91 0.2 0.1 0.05], 'visible', 'on', 'Callback', @clear_selection_func);

        addlistener(guiapp.sld_timepoint, 'Value', 'PostSet', @sliderMoving_timepoint);
        addlistener(guiapp.sld_slice, 'Value', 'PostSet', @sliderMoving_slice);

        % display
        guiapp.txt_display = uicontrol('Style', 'text', 'String', 'z: 1; t: 1', 'Units', 'normalized', 'Position', [0.005, 0, 0.15, 0.05], 'HorizontalAlignment', 'left');
        guiapp.img_display = axes('Units', 'normalized', 'Position', [0, 0.25, 1, 0.75]);
    end
    
    % general behaviour
    function keypress_callback(~, event)
        if strcmp(event.Key, 'control')
            set(f, 'WindowScrollWheelFcn', @mouse_scroll_callback_z);
        end
    end
    function keyrelease_callback(~, event)
        if strcmp(event.Key, 'control')
            set(f, 'WindowScrollWheelFcn', @mouse_scroll_callback_t);
        end
    end
    function mouse_scroll_callback_z(~, event)
        slice = slice - event.VerticalScrollCount;
        slice = max(min(slice, size(data,3)), 1);
        guiapp.sld_slice.Value = slice;
        guiapp.txt_slice.String = sprintf('Slice: %d', slice);
        update_image;
    end
    function mouse_scroll_callback_t(~, event)
        timepoint = timepoint - 3.*event.VerticalScrollCount;
        timepoint = max(min(timepoint, size(data,4)), 1);
        guiapp.sld_timepoint.Value = timepoint;
        guiapp.txt_timepoint.String = sprintf('Timepoint: %d', timepoint);
        update_image;
    end
    function resize_text_function(~,~)
        new_size=get(f,'Position');
        [~, changed_dim]=min(new_size(3:4));
        FontSize_factor=7*new_size(2+changed_dim)./current_size(2+changed_dim);
        textObjs = findobj(f, '-property', 'FontSize');
        % set(textObjs,'FontSize',FontSize_factor);
        for i=1:numel(textObjs)
            try
            textObjs(i).FontUnits="normalized";
            textObjs(i).FontSize=0.7;
            end
%             set(textObjs(i),'FontSize',FontSize_factor);
        end
        resize_spots_function();
    end
    function resize_spots_function(~,~)
        set(image_ax.Parent,'Units','pixels')
        [new_displayed_pxsz, lowest_dim]=min([image_ax.Parent.Position(3)/diff(image_ax.Parent.XLim), ...
            image_ax.Parent.Position(4)/diff(image_ax.Parent.YLim)]);
        set(image_ax.Parent,'Units','normalized');
        factor=75*(new_displayed_pxsz/(currently_displayed_pxsz(lowest_dim)))^2;
        update_image();
        %currently_displayed_pxsz=new_displayed_pxsz;
    end
    function sliderMoving_timepoint(~, ~)
        timepoint = round(guiapp.sld_timepoint.Value);
        guiapp.txt_display.String = sprintf('z: %d t: %d', slice, timepoint);
        update_image();
    end
    function sliderMoving_slice(~, ~)
        slice = round(guiapp.sld_slice.Value);
        guiapp.txt_display.String = sprintf('z: %d t: %d', slice, timepoint);
        update_image();
    end
    function color_scale_func(~, ~)
        if guiapp.bttn_caxis_control.Value==0
            caxis(image_ax.Parent,[custom_minca custom_maxca]);
            set(hRangeSlider,'Visible','on');
        else
            caxis(image_ax.Parent,'auto');
            set(hRangeSlider,'Visible','off');
        end
    end
    function hRangeSlider=caxis_slider
            caxis_control_position = get(guiapp.bttn_caxis_control,'Position');
            sliderPosition = caxis_control_position;
            sliderPosition(3)=2*sliderPosition(3); sliderPosition(4)=2*sliderPosition(4);
            sliderPosition(1) = sliderPosition(1)-sliderPosition(3);
            labelPosition = sliderPosition;
            labelPosition(2) = labelPosition(2) - labelPosition(4);
            jRS = com.jidesoft.swing.RangeSlider;
            [jRangeSlider, hRangeSlider] = javacomponent(jRS,[],f);
            set(hRangeSlider,'Units','Normalized')
            set(hRangeSlider,'Position',sliderPosition)
            set(jRangeSlider,'Maximum',maxca,...
                'Minimum',minca,...
                'LowValue',custom_minca,...
                'HighValue',custom_maxca,...
                'Name','c axis',...
                'MajorTickSpacing',round((maxca-minca)/3),...
                'MinorTickSpacing',round((maxca-minca)/30), ...
                'PaintTicks',true,...
                'PaintLabels',true,...
                'StateChangedCallback',@jRangeSlider_Callback);
            uicontrol(f,'Style','Text','Position',labelPosition,'String','c axis');
        
        function jRangeSlider_Callback(jRangeSlider,~)
            custom_minca=jRangeSlider.LowValue;
            custom_maxca=jRangeSlider.HighValue;
            custom_maxca(custom_maxca<=custom_minca)=custom_minca+1;
            caxis(image_ax.Parent,[custom_minca custom_maxca]);
        end
    end
    function clear_selection_func(~, ~)
        selection=[];
        update_plot();
    end
    % mouse pointer behaviour
    function hover_callback_removing(~, ~)
        try
            cursor_pos = get(gca,'CurrentPoint');
            x = cursor_pos(1,1);
            y = cursor_pos(1,2);
            x_data = X(:,timepoint);
            y_data = Y(:,timepoint);
            z_data = Z(:,timepoint);
            [knn_id, distance] = knnsearch([x_data, y_data],[x, y]);

            if distance < 4
                set(f,'Pointer','hand');
                hold on
                delete(hoover_over);
                delete(hoover_over_rt);
                R=125./(single(abs(z_data(knn_id)-slice)+1).^(1.5));
                hoover_over=scatter(image_ax.Parent,x_data(knn_id),y_data(knn_id),R,"red",'o','visible','on','LineWidth',3,'ButtonDownFcn',@SelectNeuron);
                hold(trace_ax,'on');
                [~,~,~,~,S]=prepareCoordData;
                hoover_over_rt=plot(trace_ax,S(knn_id,:),'Color',C(knn_id,:),'LineWidth',2);
                hold(trace_ax,'off');
                axes(image_ax.Parent);
            else
                set(gcf,'Pointer','arrow');
                delete(hoover_over);
                delete(hoover_over_rt);
            end
        end
    end
    function hover_callback_adding(~,~)
    end
    function spot_func(~,~)
        if show_spots == 1
            show_spots = 0;
        else
            show_spots = 1;
        end
        update_image();
    end

    % plot tracks
    function update_image(~,~)
        set(0, 'CurrentFigure', f);
        delete(track_plot_ax);
        delete(selected_track_plot_ax);
        delete(tm_xline);
        set(image_ax,'CData',data(:,:,slice,timepoint));
        if timepoint<=size(Z,2) && show_spots == 1
            hold on
            R=factor./(single(abs(Z(:,timepoint)-slice)+1).^(1.5));
            track_plot_ax=scatter(image_ax.Parent,X(:,timepoint),Y(:,timepoint),R,C,'LineWidth',1*sqrt(factor)/25);
            hold on
            selected_track_plot_ax=scatter(image_ax.Parent,X(selection,timepoint),Y(selection,timepoint),R(selection),'red','LineWidth',3*sqrt(factor)/25);
            tm_xline=xline(trace_ax,timepoint);
        end
        update_controls();
        drawnow
    end

    % data conversion
    function [X,Y,Z,C,S]=prepareCoordData
        X = nan(numel(current_neuron_arangement),size(current_neuron_arangement(1).coords,1));
        Y = X;
        Z = X;
        S = X;
        for i = 1:numel(current_neuron_arangement)
            for t_id=1:size(current_neuron_arangement(i).coords,1)
                X(i,t_id) = current_neuron_arangement(i).coords(t_id,2);
                Y(i,t_id) = current_neuron_arangement(i).coords(t_id,1);
                Z(i,t_id) = current_neuron_arangement(i).coords(t_id,3);
                S(i,t_id) = current_neuron_arangement(i).coords(t_id,5);
            end
        end
        cmap = colormap("jet");
        colormap("gray");
        C = cmap(round(linspace(1,256,size(X,1))),:);
    end

    % neurons selection
    function SelectNeuron(~,eventData)
        if guiapp.bttn_addneurons.Value==0
            id_x=abs(X(:,timepoint)-eventData.IntersectionPoint(1));
            id_y=abs(Y(:,timepoint)-eventData.IntersectionPoint(2));
            id_z=abs(Z(:,timepoint)-slice);
            close_points = id_x + id_y + id_z;
            [~, clicked_id] = min(close_points);
            selected_id = find(selection==clicked_id);
            if isempty(selected_id)
                selection=[selection, clicked_id];
            else
                selection(selected_id)=[];
            end
            if isempty(selection)
                set(guiapp.bttn_delneurons,'visible','off');
            else
                set(guiapp.bttn_delneurons,'visible','on');
            end
            update_image
            colormap("gray");
            axis image
            axis off

            update_plot
            disp(selection)
        end
    end
    function update_plot
        cla(trace_ax);
        if ~isempty(selection)
            [~,~,~,C,S]=prepareCoordData;
            hold(trace_ax,'on');
            for i=1:numel(selection)
                plot(trace_ax,S(selection(i),:),'Color',C(selection(i),:),'LineWidth',2);               
            end
            hold(trace_ax,'off');
        end
    end

    % remove neurons func
    function remove_neurons_func(~,~)
        to_remove=[to_remove, selection];
        update_image;
        colormap("gray");
        axis image
        axis off
    end

    % keep added location
    function adding_neuron(~,eventData)
        x=eventData.IntersectionPoint(1);
        y=eventData.IntersectionPoint(2);
        z=slice;

        to_add=cat(1,to_add,[x y z timepoint]);

        % check for tracks to track
        update_image;
        colormap("gray");
        axis image
        axis off

        disp(to_add);
    end

    % update controls when adding neurons
    function update_controls(~,~)
        if isempty(to_add)
            set(guiapp.bttn_addtracks,'visible','off');
        else
            set(guiapp.bttn_addtracks,'visible','on');
        end

        % check for changes to save
        if isempty(to_remove) && isempty(newly_tracked)
            set(guiapp.bttn_savechanges,'visible','off');
        else
            set(guiapp.bttn_savechanges,'visible','on');
        end

        % check for changes to save
        if isempty(selection)
            set(guiapp.bttn_delneurons,'visible','off');
        else
            set(guiapp.bttn_savechanges,'visible','on');
        end

        if guiapp.bttn_addneurons.Value==0
            set(f,'Pointer','arrow');
            set(image_ax,'ButtonDownFcn','');
            set(f,'WindowButtonMotionFcn', @hover_callback_removing);
        else
            set(f,'Pointer','crosshair');
            set(image_ax,'ButtonDownFcn',@adding_neuron);
            set(f,'WindowButtonMotionFcn', @hover_callback_adding);
        end
    end

    % save data
    function save_neurons_func(~,~)
        % removing neuron
        neurons_verified=current_neuron_arangement;
        neurons_verified(to_remove) = [];
        newly_tracked = [];


        % saving new neuron list
        id=previous_changes;
        save([outputpath '/neurons_reconstructed_max_' num2str(id) '.mat'],'neurons_verified','to_remove');

        % re-initializing interface
        selection=[];
        to_remove=[];
        to_add=[];
        current_neuron_arangement=neurons_verified;
        [X,Y,Z,C]=prepareCoordData;
        update_image;
        update_plot;
        colormap("gray");
        axis image
        axis off
    end

    % add_tracks_func
    function add_tracks_func(~,~)
        current_neuron_arangement=neuron_reconstruction_2(data,to_add, current_neuron_arangement,outputpath);
        [X,Y,Z,C]=prepareCoordData;
        to_add=[];
        figure(f)
        update_image;
        colormap("gray");
        axis image
        axis off
        newly_tracked=1;
    end

    % check latest file
    function latest_change=previous_changes
        previous_files=dir([outputpath '/neurons_reconstructed_max_*.mat']);
        B(numel(previous_files))=0;
        for i=1:numel(previous_files)
            B(i) = str2double(cell2mat(regexp(previous_files(i).name,'\d*','Match')));
        end
        latest_change=max(B)+1;
    end

    function stack_initialization
%         data=single(data);
        if size(data,1)>size(data,2)
            data=(permute(data,[2 1 3 4]));
        end
    end

end