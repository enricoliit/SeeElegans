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

%% creation of distmatrix3
neurons=pca_reconstruction(source_variable);
for i=1:size(added)
    for j=1:numel(neurons)
        try
            time_id=added(i,5);
            distmatrix3(i,j,1)=added(i,2)-neurons(j).coords(time_id,1);
            distmatrix3(i,j,2)=added(i,3)-neurons(j).coords(time_id,2);
            distmatrix3(i,j,3)=added(i,4)-neurons(j).coords(time_id,3);
        catch
            distmatrix3(i,j,1:3)=NaN;
        end
    end
end

%% estimated position
clear reconstruction
M=0;
for i=1:numel(neurons)
    M=max([M size(neurons(i).coords,1)]);
end

for i=1:size(added,1)
    for t_id=1:M
        total=0;
        reconstruction{i}(t_id,1:3)=[0, 0, 0];
        for j=1:numel(neurons)
            try
                if ~isnan(neurons(j).coords(t_id,1))  && ~isnan(distmatrix3(i,j,1))
                    addition=[neurons(j).coords(t_id,1)+distmatrix3(i,j,1), neurons(j).coords(t_id,2)+distmatrix3(i,j,2), neurons(j).coords(t_id,3)+distmatrix3(i,j,3)];
                    if ~isnan(addition)
                        reconstruction{i}(t_id,1:3)=reconstruction{i}(t_id,1:3)+addition;
                        total=total+1;
                    end
                end
            end
        end
        reconstruction{i}(t_id,1:3)=reconstruction{i}(t_id,1:3)./total;
    end
end

%%
a=1;
load([outputpath '/parameters.mat']);

looking_window=round(sizes/2);
looking_window_small=round(sizes/3);
radius1 =  round(2*sizes/3);
radius2= round(sizes/3+sizes/10);
radius3= round(sizes/3);
[s1, s2]=size(wholestack,[1 2]);

for add_id=1:numel(reconstruction)
    if 1
        final_time=size(reconstruction{add_id},1);
        estimated_mean_position(1:final_time,1:3)=reconstruction{add_id}(:,1:3);
        plane=round(estimated_mean_position(added(add_id,4),3));

        trace_reconstruction=NaN.*ones(1,M);
        trace_reconstruction_max=NaN.*ones(1,M);

        row=round(estimated_mean_position(1,1));
        col=round(estimated_mean_position(1,2));

        delta_t=1;

        for i=1:size(reconstruction{add_id},1)
            row=round(estimated_mean_position(i,1));
            col=round(estimated_mean_position(i,2));
            if ~isnan(row)
                row1=row-looking_window_small; row2=row+looking_window_small; col1=col-looking_window_small; col2=col+looking_window_small;
                row1(row1<1)=1; row2(row2<1)=1; col1(col1<1)=1; col2(col2<2)=1;
                row1(row1>s1)=s1; row2(row2>s1)=s1; col1(col1>s2)=s2; col2(col2>s2)=s2;

                filtered_image=imgaussfilt(wholestack(row1:row2,col1:col2,plane,i),1);
                [row_max,col_max]=find(filtered_image==max(max(filtered_image)));
                row_max=row_max(1); col_max=col_max(1);

                row_1i=row-looking_window_small; row_2i=row+looking_window_small; col_1i=col-looking_window_small; col_2i=col+looking_window_small;
                row_1i(row_1i<1)=1; col_1i(col_1i<1)=1;
                row_2i(row_2i>s1)=s1; col_2i(col_2i>s2)=s2;

                row_max1=row_1i; row_max2=row_2i;
                col_max1=col_1i; col_max2=col_2i;
                row_max1(row_max1<1)=1; row_max2(row_max2<1)=1; col_max1(col_max1<1)=1; col_max2(col_max2<1)=1;
                row_max1(row_max1>s1)=s1; row_max2(row_max2>s1)=s1; col_max2(col_max2>s2)=s2; col_max1(col_max1>s2)=s2;

                background=mean(mean(wholestack(8:36,8:36,plane,i)'));
                trace_reconstruction(i)=mean(mean((wholestack(row1:row2,col1:col2,plane,i)-background)./1));
                %
                row_c=row1+row_max-1; col_c=col1+col_max-1;

                % CORONA ----------------------------------------------
                [columnsInImage, rowsInImage] = meshgrid(1:s2, 1:s1);
                center_col = col_c;
                center_row = row_c;
                c_col(i,1)=center_col;
                c_row(i,1)=center_row;

                circlePixels = (rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 <= radius1.^2 & (rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 >= radius2.^2;
                circlePixels_in = ((rowsInImage - center_row).^2 + (columnsInImage - center_col).^2 <= radius3.^2);

                cP=(circlePixels(row_1i:row_2i,col_1i:col_2i));
                cP_in=(circlePixels_in(row_1i:row_2i,col_1i:col_2i));

                try
                    ROI_u=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane+1,i:i+3),3)-background)./1,1);
                catch
                    ROI_u=NaN.*ones(size(ROI_m));
                end
                ROI_m=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane,i:i+3),3)-background)./1,1);
                try
                    ROI_d=imgaussfilt((median(wholestack(row_max1:row_max2,col_max1:col_max2,plane-1,i:i+3),3)-background)./1,1);
                catch
                    ROI_d=NaN.*ones(size(ROI_m));
                end
                pixels_from_ROIs=[ROI_u(cP_in); ROI_m(cP_in); ROI_d(cP_in)];
                chosen_pixels=sort(pixels_from_ROIs,'descend');
                threshold_on_pixels=round(numel(chosen_pixels)*0.1)+1;

                chosen_pixels(threshold_on_pixels:end)=[];
                chosen_pixels(isnan(chosen_pixels))=[];

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

        title([num2str(add_id)]);
        trace_reconstruction(isnan(trace_reconstruction))=[];
        timepoints=numel(trace_reconstruction);
        neurons_reconstructed_added(add_id).coords(1:timepoints,5)=trace_reconstruction(1:timepoints);
        neurons_reconstructed_added(add_id).coords(1:timepoints,[1 2])=[(estimated_mean_position(1:timepoints,1)), (estimated_mean_position(1:timepoints,2))];
        neurons_reconstructed_added(add_id).coords(1:timepoints,[3 4])=[estimated_mean_position(1:timepoints,3), (1:timepoints)'];
        neurons_reconstructed_max_added(add_id).coords(1:timepoints,5)=(trace_reconstruction_max(1:timepoints));
        neurons_reconstructed_max_added(add_id).coords(1:timepoints,[1 2])=[(c_row(1:timepoints,:)), (c_col(1:timepoints,:))];
        neurons_reconstructed_max_added(add_id).coords(1:timepoints,[3 4])=[estimated_mean_position(1:timepoints,3), (1:timepoints)'];

        a=a+1;
        if a>50
            a=1;
            figure
        end
    end
end
file_to_convert=neurons_reconstructed_added;
output_variable=source_variable;

for i=1:numel(file_to_convert)
    output_variable(end+1).coords=file_to_convert(i).coords;
end

end


