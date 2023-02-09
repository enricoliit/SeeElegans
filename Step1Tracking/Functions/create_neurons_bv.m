function neurons_bv=create_neurons_bv(short_neurons_tmp_unsorted_s,reconstructed_links)
new_traces{1}=[];
if ~isempty(reconstructed_links)
    new_traces{numel(reconstructed_links)}=[];
    for i=1:numel(reconstructed_links)
        new_traces{i}.coords=[];
        for j=1:numel(reconstructed_links{i})
            new_traces{i}.coords=[new_traces{i}.coords; short_neurons_tmp_unsorted_s(reconstructed_links{i}(j)).coords];
        end
    end
    % merging linked and complete traces
    skip=unique([reconstructed_links{:}]);
    a=1;
    neurons_bv=[];
    for i=1:size(short_neurons_tmp_unsorted_s,2)
        if ismember(i,skip)==0
            neurons_bv(a).coords=short_neurons_tmp_unsorted_s(i).coords;
            a=a+1;
        end
    end
    for i=1:numel(new_traces)
        neurons_bv(a).coords=new_traces{i}.coords;
        a=a+1;
    end
else
    a=1;
    neurons_bv=[];
    for i=1:size(neurons_tmp_unsorted_s,2)
        neurons_bv(a).coords=short_neurons_tmp_unsorted_s(i).coords;
        a=a+1;
    end
end
