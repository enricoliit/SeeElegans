function output_variable=neuron_reconstruction(wholestack,to_be_added,source_variable,outputpath)
% wholestack: the recorded stacked
% to_be_added: the neurons to be added
% source_variable: neurons already tracked
% outputpath: outputfolder

figure
id_added(size(to_be_added,1),1)=0;
for i=1:size(to_be_added,1)
    id_added(i,1)=i;
end
added=cat(2,id_added,to_be_added);
added=added(:,[1 3 2 4 5]);
neurons=pca_reconstruction(source_variable);

% creation of distmatrix3
% distMatrix3=makeDistMatrix(added,neurons);

clear reconstruction
t_max = max(cellfun(@(x) size(x, 1), {neurons.coords}));

for add_id = 1:size(added,1)

    for j = 1:numel(neurons)
        distMatrix3(add_id,j,:)=distMatrixElmnt(added(add_id,:),neurons(j));

        try
            velMatrix3(add_id,j,1)=added(add_id,2)-neurons(j).coords(time_id,1);
            velMatrix3(add_id,j,2)=added(add_id,3)-neurons(j).coords(time_id,2);
            velMatrix3(add_id,j,3)=added(add_id,4)-neurons(j).coords(time_id,3);
        catch
            velMatrix3(add_id,j,1:3)=NaN;
        end
    end

    for t_id=1:t_max
        reconstruction{add_id}(t_id,1:3)=reconstructTrack(neurons, distMatrix3, add_id, t_id);
    end
end

%%
a=1;
load([outputpath '/parameters.mat']);
search_factor = 0.5;                                                       % 0.5 - 3
looking_window = round(search_factor .* 0.85 .* sizes);                    % SCALA
looking_window_small = round(search_factor .* 0.6 .* sizes);               % SCALA

% parameters to measure intensity
radius1 = sizes;                                                           % SCALA
radius2 = round(0.7*sizes);                                                % SCALA
% radius3 = round(0.6*sizes);                                              % SCALA
radius3 = looking_window_small;

[s1, s2]=size(wholestack,[1 2]);

for add_id=1:size(added,1)
    if 1
        final_time=size(reconstruction{add_id},1);
        estimated_mean_position(1:final_time,1:3)=reconstruction{add_id}(:,1:3);
        plane=round(estimated_mean_position(added(add_id,4),3));

        trace_reconstruction=NaN.*ones(1,t_max);
        trace_reconstruction_max=NaN.*ones(1,t_max);

        row=round(estimated_mean_position(1,1));
        col=round(estimated_mean_position(1,2));

        delta_t=1;
        iRoi = 1;
        for t_id=1:delta_t:size(reconstruction{add_id},1)
            
            row = round(estimated_mean_position(t_id,1));
            col = round(estimated_mean_position(t_id,2));
            
            if ~isnan(row)

                % looking window small extrema
                row1=row-looking_window_small; row2=row+looking_window_small; col1=col-looking_window_small; col2=col+looking_window_small;
                row1(row1<1)=1; row2(row2<1)=1; col1(col1<1)=1; col2(col2<2)=1;
                row1(row1>s1)=s1; row2(row2>s1)=s1; col1(col1>s2)=s2; col2(col2>s2)=s2;

                % finding maximum value in looking window small
                filtered_image=imgaussfilt(wholestack(row1:row2,col1:col2,plane,t_id),1); % SCALA
                [row_max,col_max]=find(filtered_image==max(max(filtered_image)));
                row_max=row_max(1); col_max=col_max(1);

                % looking window extrema
                row_1i=row-looking_window; row_2i=row+looking_window; col_1i=col-looking_window; col_2i=col+looking_window;
                row_1i(row_1i<1)=1; col_1i(col_1i<1)=1;
                row_2i(row_2i>s1)=s1; col_2i(col_2i>s2)=s2;

                % duplicate extrema
                row_max1=row_1i; row_max2=row_2i;
                col_max1=col_1i; col_max2=col_2i;
                row_max1(row_max1<1)=1; row_max2(row_max2<1)=1; col_max1(col_max1<1)=1; col_max2(col_max2<1)=1;
                row_max1(row_max1>s1)=s1; row_max2(row_max2>s1)=s1; col_max2(col_max2>s2)=s2; col_max1(col_max1>s2)=s2;

                % raw intensity background subtracted
                background=mean(mean(wholestack(8:36,8:36,plane,t_id)'));
                trace_reconstruction(t_id)=mean(mean((wholestack(row1:row2,col1:col2,plane,t_id)-background)./1));

                % row and col of max intensity in whole slice
                row_c=row1+row_max-1; col_c=col1+col_max-1;

                % CORONA ----------------------------------------------
                [columnsInImage, rowsInImage] = meshgrid(1:s2, 1:s1);
                center_col = col_c;
                center_row = row_c;
                c_col(t_id,1) = center_col;
                c_row(t_id,1) = center_row;

                circlePixels = (rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 <= radius1.^2 & (rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 >= radius2.^2;
                circlePixels_in = ((rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 <= radius3.^2);

                cP = (circlePixels(row_1i:row_2i,col_1i:col_2i));
                cP_in = (circlePixels_in(row_1i:row_2i,col_1i:col_2i));
                ROI_m = imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane,t_id:t_id+3),4)-background)./1,1);
                %%-------------------------------------------------------------------------
                ROIs{iRoi}=[ROI_m];
                iRoi=iRoi+1;
                %%-------------------------------------------------------------------------

                try
                    ROI_u = imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane+1,t_id:t_id+3),4)-background)./1,1);
                catch
                    ROI_u = NaN.*ones(size(ROI_m));
                end
                try
                    ROI_d = imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane-1,t_id:t_id+3),4)-background)./1,1);
                catch
                    ROI_d = NaN.*ones(size(ROI_m));
                end
                pixels_from_ROIs = [ROI_u(cP_in); ROI_m(cP_in); ROI_d(cP_in)];
                chosen_pixels = sort(pixels_from_ROIs,'descend');
                threshold_on_pixels = round(numel(chosen_pixels)*0.1)+1;

                chosen_pixels(threshold_on_pixels:end) = [];
                chosen_pixels(isnan(chosen_pixels)) = [];

                if ~isempty((chosen_pixels))
                    if ~isempty(ROI_m(find(cP)))
                        trace_reconstruction_max(t_id) = median(chosen_pixels)-median(ROI_m(find(cP)));
                    else
                        trace_reconstruction_max(t_id) = NaN;
                    end
                    trace_reconstruction_max_no_corona(t_id) = median(chosen_pixels);
                else
                    trace_reconstruction_max_no_corona(t_id) = NaN;
                end
            else
                trace_reconstruction_max(t_id) = NaN;
                trace_reconstruction_max_no_corona(t_id) = NaN;
            end
        end

        
        subplot(5,10,a)
        set(gca, 'ColorOrderIndex', 1) %
        trel=numel(trace_reconstruction(1:delta_t:end));
        p1=plot(1:delta_t:trel,(trace_reconstruction(1:delta_t:trel)));
        hold on
        p12=plot(1:delta_t:t_id,(trace_reconstruction_max(1:delta_t:t_id)));
        hold on
        p15=plot(1:delta_t:t_id,(trace_reconstruction_max_no_corona(1:delta_t:t_id)));

        title([num2str(add_id)]);
        trace_reconstruction(isnan(trace_reconstruction))=[];
        timepoints=numel(trace_reconstruction);
        neurons_reconstructed_added(add_id).coords(1:timepoints,5)=trace_reconstruction(1:timepoints);
        neurons_reconstructed_added(add_id).coords(1:timepoints,[1 2])=[(estimated_mean_position(1:timepoints,1)), (estimated_mean_position(1:timepoints,2))];
        neurons_reconstructed_added(add_id).coords(1:timepoints,[3 4])=[estimated_mean_position(1:timepoints,3), (1:timepoints)'];
        % ATTENTION
        neurons_reconstructed_max_added(add_id).coords(1:timepoints,5)=(trace_reconstruction(1:timepoints));
        % ---------
        neurons_reconstructed_max_added(add_id).coords(1:timepoints,[1 2])=[(c_row(1:timepoints,:)), (c_col(1:timepoints,:))];
        neurons_reconstructed_max_added(add_id).coords(1:timepoints,[3 4])=[estimated_mean_position(1:timepoints,3), (1:timepoints)'];

        a=a+1;
        if a>50
            a=1;
            figure
        end
    end
end
file_to_convert=neurons_reconstructed_max_added;
output_variable=source_variable;

for i=1:numel(file_to_convert)
    output_variable(end+1).coords=file_to_convert(i).coords;
end



end

function distMatrix3=makeDistMatrix(added,neurons)
distMatrix3 = zeros(size(added,1),numel(neurons),3);
for i=1:size(added,1)
    for j=1:numel(neurons)
        distMatrix3(i,j,:) = singleDistMatrix(i,j,added,neurons);
    end
end
end

function distMatrix=singleDistMatrix(i,j,added,neurons)
% estimated position
try
    time_id=added(i,5);
    distMatrix(1)=added(i,2)-neurons(j).coords(time_id,1);
    distMatrix(2)=added(i,3)-neurons(j).coords(time_id,2);
    distMatrix(3)=added(i,4)-neurons(j).coords(time_id,3);
catch
    distMatrix(1:3)=NaN;
end
end

function distMatrix=distMatrixElmnt(added,neurons)
% estimated position
try
    time_id=added(1,5);
    distMatrix(:)=added(1,2:4)-neurons.coords(time_id,1:3);
%     distMatrix(1)=added(1,2)-neurons.coords(time_id,1);
%     distMatrix(2)=added(1,3)-neurons.coords(time_id,2);
%     distMatrix(3)=added(1,4)-neurons.coords(time_id,3);
catch
    distMatrix(1:3)=NaN;
end
end

function reconstruction = reconstructTrack(neurons, distMatrix3, i, t_id)

%         total = 0;
        total_weights = 0;
        reconstruction=[0, 0, 0];
        for j=1:numel(neurons)
            if size(neurons(j).coords,1) >= t_id 
                if ~isnan(neurons(j).coords(t_id,1))  && ~isnan(distMatrix3(i,j,1))
                    addition=neurons(j).coords(t_id,1:3) + squeeze(distMatrix3(i,j,:))';
                    single_weight = sum(distMatrix3(i,j,:).^2);
                    if ~isnan(addition)
%                         reconstruction{i}(t_id,1:3) = reconstruction{i}(t_id,1:3) + addition;
                        reconstruction = reconstruction + addition .* single_weight; % weighted mean
                        total_weights = total_weights + single_weight;
%                         total = total+1;
                    end
                end
            end
        end
%         reconstruction{i}(t_id,1:3)=reconstruction{i}(t_id,1:3)./total;
        reconstruction = reconstruction./total_weights;
end