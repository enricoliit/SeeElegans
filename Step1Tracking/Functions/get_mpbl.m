function mutual_position_based_links=get_mpbl(distmatrix_s,short_neurons_tmp_unsorted_s,maxtime_s,xy_pixel_size,z_pixel_size,max_gap_closing_mutual_distance)

distmatrix_s_mean=squeeze(mean(distmatrix_s,3,'omitnan'));
taken=[];
global p, global h; D = parallel.pool.DataQueue; h = waitbar(0, 'Merging tracks based on mutual distance ...'); afterEach(D, @nUpdateWaitbar_s); p=1;
for i=1:numel(short_neurons_tmp_unsorted_s)
    send(D,i)
    if short_neurons_tmp_unsorted_s(i).coords(end,4)<maxtime_s
        link=0;
        time_constraint=short_neurons_tmp_unsorted_s(i).coords(end,4);
        loop_set=1:size(distmatrix_s,1); loop_set(loop_set==i)=[]; loop_set(ismember(loop_set,taken))=[];
        best_candidate=Inf;
        saved=NaN*ones(1,max(loop_set));
        for j=loop_set
            if short_neurons_tmp_unsorted_s(j).coords(1,4)>=time_constraint
                current_difference=sqrt( xy_pixel_size.*(distmatrix_s_mean(i,:,1)-distmatrix_s_mean(j,:,1)).^2 + xy_pixel_size.*(distmatrix_s_mean(i,:,2)-distmatrix_s_mean(j,:,2)).^2 + z_pixel_size.*(distmatrix_s_mean(i,:,3)-distmatrix_s_mean(j,:,3)).^2); % aggiungo il tempo?
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
close(h)
