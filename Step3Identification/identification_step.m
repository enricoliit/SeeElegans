function identification_step(data,neurons_cleaned,outputpath,voxel_size)
% 

addpath('./Functions/');

% variables initialization
timepoint = 1;
slice = 1;
selection=[];
assigning=0;

% directions=[];
chosen_point=[];
clicked_id_to_assign=[];
if nargin<4
    voxel_size=[0.267 0.267 2];
end
zratio=voxel_size(3)/voxel_size(1);

% neuron names 
namesstr={'RMEL','RMER','RMEV','RMED','RID','RIS','RIBL','RIBR','AVBL',...
    'AVBR','VB1','VB2','VB3','DB1','DB2','SIAVL','SIAVR','AVAL','AVAR',...
    'AVEL','AVER','AIBL','AIBR','RIML','RIMR','SMDVL','SMDVR','RIVL',...
    'RIVR','SMDDL','SMDDR','RMDL','RMDR','RMDDL','RMDDR','RMDVL'};
current_neuron_arangement = pca_reconstruction(neurons_cleaned,20);
current_id_hypothesis=[];
current_id_hypothesis_string{1}=[];
manual_id=[];
[X,Y,Z,C]=prepareCoordData;

% controls
f=figure('Units','normalized','visible','off');
guiapp.sld_timepoint = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(data,4), 'Value', 1, 'Units', 'normalized', 'Position', [0.16, 0, 0.72, 0.04], 'Callback', @sliderMoving_timepoint);
guiapp.sld_slice = uicontrol('Style', 'slider', 'Min', 1, 'Max', size(data,3), 'SliderStep', [1/(size(data,3)-1) 4/(size(data,3)-1)], 'Value', 1, 'Units', 'normalized', 'Position', [0.16, 0.04, 0.72, 0.04], 'Callback', @sliderMoving_slice);
guiapp.bttn_manualassign = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Assign', 'FontSize', 6, 'Position', [0.005, 0.05, 0.1, 0.05], 'visible', 'on', 'Callback', @manual_identity_func);
guiapp.bttn_autoassign = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Auto Assign', 'FontSize', 6, 'Position', [0.005, 0.1, 0.1, 0.05], 'visible', 'on', 'Callback', @automatic_identity_func);
guiapp.bttn_headpoint = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Head Point', 'FontSize', 6, 'Position', [0.005, 0.15, 0.1, 0.05], 'visible', 'on', 'Callback', @choose_head_point);
guiapp.text_selassign = uicontrol('Style', 'text', 'Units', 'normalized', 'String', 'id: ns', 'FontSize', 6, 'Position', [0.895, 0.23, 0.1, 0.05], 'visible', 'on', 'Callback', @choose_head_point);
guiapp.popup_idoptions = uicontrol('Style', 'popupmenu', 'Units', 'normalized', 'String', namesstr, 'FontSize', 6, 'Position', [0.895, 0.135, 0.1, 0.05], 'visible', 'on');
guiapp.bttn_storeid = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Store', 'FontSize', 6, 'Position', [0.92, 0.07, 0.08, 0.05], 'visible', 'on', 'Callback', @modify_current_hypothesis);
guiapp.bttn_unassign = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Delete', 'FontSize', 6, 'Position', [0.92, 0.2, 0.08, 0.05], 'visible', 'on', 'Callback', @unassign);
guiapp.bttn_savechanges = uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'String', 'Save', 'FontSize', 6, 'Position', [0.92, 0.01, 0.08, 0.05], 'visible', 'on', 'Callback', @save_neurons_func);
addlistener(guiapp.sld_timepoint, 'Value', 'PostSet', @sliderMoving_timepoint);
addlistener(guiapp.sld_slice, 'Value', 'PostSet', @sliderMoving_slice);

% display axes
guiapp.txt_display = uicontrol('Style', 'text', 'String', 'z: 1; t: 1', 'Units', 'normalized', 'Position', [0.005, 0, 0.15, 0.05], 'HorizontalAlignment', 'left');
guiapp.img_display = axes('Units', 'normalized', 'Position', [0, 0.25, 1, 0.75/2]);
guiapp.s3D_display = axes('Units', 'normalized', 'Position', [0, 0.25+0.75/2, 1, 0.75/2]);
guiapp.trc_display = axes('Units','normalized','Position',[0.2 0.12 0.65 0.11]);
hold(guiapp.img_display,'all');
hold(guiapp.s3D_display,'all');
hold(guiapp.trc_display,'all');

% GUI elements initialization
hoover_over=plot(0,0,'.','Visible','off');
text_on_3D_plot=text(guiapp.img_display,1,1,1,'');
text_on_slice=text(guiapp.s3D_display,1,1,1,'');

% update the image display
track_plot_ax=scatter(guiapp.img_display,X(:,timepoint),Y(:,timepoint),1,'or');
hold(guiapp.img_display,'all'); % da levare
selected_track_plot_ax=scatter(guiapp.img_display,NaN,NaN);
image_plot=imagesc(guiapp.img_display,data(:,:,slice,timepoint));
s3D_plot=scatter3(guiapp.s3D_display,X(:,timepoint),Y(:,timepoint),Z(:,timepoint).*zratio,1,'red');
hold on
hs_lim(1,:)=get(guiapp.s3D_display,'Xlim');
hs_lim(2,:)=get(guiapp.s3D_display,'Ylim');
hs_edges_x=[hs_lim(1,1) hs_lim(1,2) hs_lim(1,2) hs_lim(1,1)];
hs_edges_y=[hs_lim(2,1) hs_lim(2,1) hs_lim(2,2) hs_lim(2,2)];
hs_edges_z=[slice slice slice slice].*zratio;
slice_highlight_plot=fill3(guiapp.s3D_display,hs_edges_x, hs_edges_y, hs_edges_z,'red','FaceAlpha',0.5);
hold on
colormap gray;
axis(guiapp.s3D_display,'equal','off');
grid(guiapp.s3D_display,'off');
xlim(guiapp.s3D_display,'manual');
axis(guiapp.img_display,'equal','off');
axis(guiapp.img_display);
tm_xline=xline(guiapp.trc_display,timepoint);
axes(guiapp.img_display);
set(f, 'Visible', 'on');
update_image

% controls
set(f, 'WindowButtonMotionFcn', @hover_callback_showing);
set(f, 'WindowScrollWheelFcn', @mouse_scroll_callback_t);
set(f, 'Windowkeypressfcn', @keypress_callback);
set(f, 'WindowKeyReleaseFcn',@keyrelease_callback);

% ------------------------- functions -------------------------------------

% mouse controls ------------------------------------------------------
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
        timepoint = timepoint - event.VerticalScrollCount;
        timepoint = max(min(timepoint, size(data,4)), 1);
        guiapp.sld_timepoint.Value = timepoint;
        guiapp.txt_timepoint.String = sprintf('Timepoint: %d', timepoint);
        update_image;
    end
    function hover_callback_showing(~, ~)
        if timepoint<=size(Z,2)
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
                R=125./(single(abs(z_data(knn_id)-slice)+1).^(1.5));
                hoover_over=scatter(guiapp.img_display,x_data(knn_id),y_data(knn_id),R,"red",'o','visible','on','LineWidth',3,'ButtonDownFcn',@SelectNeuron);
                axes(guiapp.img_display);
            else
                set(gcf,'Pointer','arrow');
                delete(hoover_over);
            end
        else
            set(gcf,'Pointer','arrow');
            delete(hoover_over);
        end
    end
    function hover_callback_assigning(~, ~)
        if timepoint<=size(Z,2)
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
                R=125./(single(abs(z_data(knn_id)-slice)+1).^(1.5));
                hoover_over=scatter(guiapp.img_display,x_data(knn_id),y_data(knn_id),R,"yellow",'o','visible','on','LineWidth',3,'ButtonDownFcn',@being_assigned);
            else
                set(gcf,'Pointer','arrow');
                delete(hoover_over);
            end
        else
            set(gcf,'Pointer','arrow');
            delete(hoover_over);
        end
    end

% button controls -----------------------------------------------------
    function sliderMoving_timepoint(~, ~)
        timepoint = round(guiapp.sld_timepoint.Value);
        guiapp.txt_display.String = sprintf('z: %d t: %d', slice, timepoint);
        update_image
    end
    function sliderMoving_slice(~, ~)
        slice = round(guiapp.sld_slice.Value);
        guiapp.txt_display.String = sprintf('z: %d t: %d', slice, timepoint);
        update_image
    end

% computation for graphics --------------------------------------------
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
        if  assigning==0
            cmap = colormap("jet");
            colormap("gray");
            C = cmap(round(linspace(1,256,size(X,1))),:);
        else
            cmap = colormap("jet");
            colormap("gray");
            corrM = correlation_matrix(current_neuron_arangement);
            corr_clckd_id=corrM(clicked_id_to_assign,:);
            col_corr_order=round(mat2gray(corr_clckd_id).*255)+1;
            C = cmap(col_corr_order,:);
        end
    end

% GUI elements update -------------------------------------------------
    function update_image(~,~)
        set(0, 'CurrentFigure', f);
        delete(track_plot_ax);
        delete(selected_track_plot_ax);

        % update slice plot
        set(image_plot,'CData',data(:,:,slice,timepoint));
        if timepoint<=size(Z,2)
            % update ROIs and text in slice plot
            R=125./(single(abs(Z(:,timepoint)-slice)+1).^(1.5));
            track_plot_ax=scatter(guiapp.img_display,X(:,timepoint),Y(:,timepoint),R,C,'LineWidth',1);
            hold on
            selected_track_plot_ax=scatter(guiapp.img_display,X(selection,timepoint),Y(selection,timepoint),R(selection),C(selection,:),'LineWidth',3);
            hold on
            selected_track_plot_ax(end+1)=scatter(guiapp.img_display,X(clicked_id_to_assign,timepoint),Y(clicked_id_to_assign,timepoint),R(clicked_id_to_assign),'yellow','LineWidth',3);
            hold on
            % update ROIs and text in 3D reconstruction
            set(s3D_plot,'XData',X(:,timepoint),'YData',Y(:,timepoint),'ZData',Z(:,timepoint).*zratio,'SizeData',5,'CData',C);
            if ~isempty(current_id_hypothesis)
                update_text_3D
            end
        end
        highlight_plane
        set(tm_xline,'Value',timepoint);
        update_controls
        axes(guiapp.img_display);
        drawnow
    end

    function update_text_3D
        for i=1:numel(current_id_hypothesis)
            %             text_on_slice(i)=text(guiapp.s3D_display,X(current_id_hypothesis(i),timepoint), Y(current_id_hypothesis(i),timepoint),'String',current_id_hypothesis_string{i});
            set(text_on_3D_plot(i),'Position',[X(current_id_hypothesis(i),timepoint), Y(current_id_hypothesis(i),timepoint), Z(current_id_hypothesis(i),timepoint).*zratio],'String',current_id_hypothesis_string{i});
        end

        for i=1:numel(current_id_hypothesis)
            if abs(Z(current_id_hypothesis(i),timepoint)-slice)<=1.5
                set(text_on_slice(i),'Position',[X(current_id_hypothesis(i),timepoint), Y(current_id_hypothesis(i),timepoint)],'String',current_id_hypothesis_string{i},'visible','on','color',[1 1 1],'FontSize',8);
            else
                set(text_on_slice(i),'visible','off');
            end
        end

    end

    function update_controls(~,~)
        if ~guiapp.bttn_headpoint.Value
            set(f,'Pointer','arrow');
            set(f,'WindowButtonMotionFcn', @hover_callback_showing);
        else
            choose_head_point
        end
    end

    function choose_head_point(~,~)
        set(f,'WindowButtonMotionFcn','');
        set(f,'Pointer','crosshair');
        set(image_plot,'ButtonDownFcn',@chosen_head_point);
    end

    function SelectNeuron(~,eventData)
        if guiapp.bttn_autoassign.Value==0
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
            update_image
            colormap("gray");
            axis image
            axis off

            update_plot
            axes(guiapp.img_display);
            disp(selection)
        end
    end

    function update_plot
        cla(guiapp.trc_display);
        tm_xline=xline(guiapp.trc_display,timepoint);
        if ~isempty(selection)
            [~,~,~,C,S]=prepareCoordData;
            hold(guiapp.trc_display,'on');
            for i=1:numel(selection)
                plot(guiapp.trc_display,mat2gray(S(selection(i),:)),'Color',C(selection(i),:));
            end
            hold(guiapp.trc_display,'off');
        end
    end

    function highlight_plane
        set(slice_highlight_plot,'ZData',[slice slice slice slice].*zratio);
    end

    function initialize_id_texts(~,~)
        delete(text_on_slice);
        delete(text_on_3D_plot);
        for i=1:numel(current_id_hypothesis)
            text_on_3D_plot(i)=text(guiapp.s3D_display,1,1,1,'');
        end
        for i=1:numel(current_id_hypothesis)
            text_on_slice(i)=text(guiapp.img_display,1,1,1,'');
        end
    end

    function automatic_identity_func(~,~)
        [correspondences,directions]=automatic_identity(current_neuron_arangement,chosen_point,voxel_size,round(size(data,4)/2));
        delete(text_on_slice);
        delete(text_on_3D_plot);
        if ~isempty(correspondences)
            current_id_hypothesis=[correspondences{:,1}]';
            current_id_hypothesis_string=correspondences(:,2);
            disp(current_id_hypothesis)
            initialize_id_texts;
        end
        update_image;
    end

    function manual_identity_func(~,~)
        assigning=1;
        clicked_id_to_assign=[];
        set(f, 'WindowButtonMotionFcn', @hover_callback_assigning);
    end

    function being_assigned(~,eventData)
        id_x=abs(X(:,timepoint)-eventData.IntersectionPoint(1));
        id_y=abs(Y(:,timepoint)-eventData.IntersectionPoint(2));
        id_z=abs(Z(:,timepoint)-slice);
        close_points = id_x + id_y + id_z;
        [~, clicked_id_to_assign] = min(close_points);
        guiapp.text_selassign.String = sprintf('id: %d', clicked_id_to_assign);
        [~,~,~,C]=prepareCoordData;
        update_image;
        
        set(f, 'WindowButtonMotionFcn', @hover_callback_showing);
    end

    function chosen_head_point(~,eventData)
        chosen_point=eventData.IntersectionPoint;
        update_controls;
        axes(guiapp.img_display);
    end

    function modify_current_hypothesis(~,~)
        if ~isempty(clicked_id_to_assign)
            to_be_updated=find(current_id_hypothesis==clicked_id_to_assign);
            manual_id(manual_id==clicked_id_to_assign)=[];
            current_id_hypothesis(to_be_updated)=[];
            current_id_hypothesis_string(to_be_updated)=[];
            current_id_hypothesis(end+1)=clicked_id_to_assign;
            manual_id(end+1)=clicked_id_to_assign;
            if ~isempty(current_id_hypothesis_string)
                if isempty(current_id_hypothesis_string{1})
                    current_id_hypothesis_string{1}=namesstr(guiapp.popup_idoptions.Value);
                else
                    current_id_hypothesis_string{end+1}=namesstr(guiapp.popup_idoptions.Value);
                end
            else
                current_id_hypothesis_string{1}=namesstr(guiapp.popup_idoptions.Value);
            end
            clicked_id_to_assign=[];
            assigning=0;
            [~,~,~,C,~]=prepareCoordData;
            initialize_id_texts
            update_image;
        end
    end

    function unassign(~,~)
        if ~isempty(clicked_id_to_assign)
            to_be_updated=find(current_id_hypothesis==clicked_id_to_assign);
            manual_id(manual_id==clicked_id_to_assign)=[];
            current_id_hypothesis(to_be_updated)=[];
            current_id_hypothesis_string(to_be_updated)=[];

            clicked_id_to_assign=[];
            initialize_id_texts
            update_image;
        end
    end

    % save data
    function save_neurons_func(~,~)
        % removing neuron
        neurons_identified.id=current_id_hypothesis;
        neurons_identified.name=current_id_hypothesis_string;
        % saving new neuron list
        id=previous_id_changes;
        save([outputpath '\neurons_reconstructed_max_id_' num2str(id) '.mat'],'neurons_cleaned','neurons_identified');
    end

    % check latest file
    function latest_change=previous_id_changes
        previous_files=dir('.\neurons_reconstructed_max_*.mat');
        if numel(previous_files)>0
            B(numel(previous_files))=0;
            for i=1:numel(previous_files)
                B(i) = str2double(cell2mat(regexp(previous_files(i).name,'\d*','Match')));
            end
        else
            B(1)=0;
        end
        latest_change=max(B)+1;
    end
end