function distmatrix_s=makeDistMatrix_s(neurons_tmp_unsorted_s,maxtime_s)
h=waitbar(0,'Calculating distance matrix...');
distmatrix_s=NaN*ones(numel(neurons_tmp_unsorted_s),numel(neurons_tmp_unsorted_s),1,3);
for i=1:numel(neurons_tmp_unsorted_s)
    waitbar(i/numel(neurons_tmp_unsorted_s),h)
    valid_times=neurons_tmp_unsorted_s(i).coords(:,4);
    logic_times=~isnan(valid_times);
    valid_times_ids=valid_times(logic_times);
    for k=1:numel(neurons_tmp_unsorted_s)
        valid_times_k=neurons_tmp_unsorted_s(k).coords(:,4);
        logic_times_k=~isnan(valid_times_k);
        valid_times_ids_k=valid_times_k(logic_times_k);
        i_matrix=NaN.*ones(maxtime_s,3); i_matrix(valid_times_ids,:)=neurons_tmp_unsorted_s(i).coords(logic_times,1:3);
        k_matrix=NaN.*ones(maxtime_s,3); k_matrix(valid_times_ids_k,:)=neurons_tmp_unsorted_s(k).coords(logic_times_k,1:3);
        distmatrix_s(i,k,1,1:3)= mean((i_matrix - k_matrix),1,'omitnan');
    end
end
close(h)