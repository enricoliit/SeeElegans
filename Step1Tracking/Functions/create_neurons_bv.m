function neurons_bv=create_neurons_bv(coords_t,reconstructed_links)
new_tracks{1}=[];
if ~isempty(reconstructed_links)
    new_tracks{numel(reconstructed_links)}=[];
    for i=1:numel(reconstructed_links)
        new_tracks{i}.coords=[];
        for j=1:numel(reconstructed_links{i})
            segment=squeeze(coords_t(:,:,reconstructed_links{i}(j)));
            segment = segment(~any(isnan(segment), 2), :);
            new_tracks{i}.coords=[new_tracks{i}.coords; segment];
        end
    end
    % merging linked and complete traces
    skip=unique([reconstructed_links{:}]);
    
    a=1;
    neurons_bv(size(coords_t,3)+numel(new_tracks)-numel(skip)).coords=[];
    for i=1:size(coords_t,3)
        if ismember(i,skip)==0
            segment=squeeze(coords_t(:,:,i));
            segment = segment(~any(isnan(segment), 2), :);
            neurons_bv(a).coords=segment;
            a=a+1;
        end
    end
    for i=1:numel(new_tracks)
        neurons_bv(a).coords=new_tracks{i}.coords;
        a=a+1;
    end
else
    a=1;
    neurons_bv(size(coords_t,3))=[];
    for i=1:size(coords_t,2)
        segment=squeeze(coords_t(:,:,i));
        segment = segment(~any(isnan(segment), 2), :);
        neurons_bv(a).coords=segment;
        a=a+1;
    end
end


% function neurons_bv=create_neurons_bv(short_neurons_tmp_unsorted_s,reconstructed_links)
% new_tracks{1}=[];
% if ~isempty(reconstructed_links)
%     new_tracks{numel(reconstructed_links)}=[];
%     for i=1:numel(reconstructed_links)
%         new_tracks{i}.coords=[];
%         for j=1:numel(reconstructed_links{i})
%             new_tracks{i}.coords=[new_tracks{i}.coords; short_neurons_tmp_unsorted_s(reconstructed_links{i}(j)).coords];
%         end
%     end
%     % merging linked and complete traces
%     skip=unique([reconstructed_links{:}]);
%     a=1;
%     neurons_bv=[];
%     for i=1:size(short_neurons_tmp_unsorted_s,2)
%         if ismember(i,skip)==0
%             neurons_bv(a).coords=short_neurons_tmp_unsorted_s(i).coords;
%             a=a+1;
%         end
%     end
%     for i=1:numel(new_tracks)
%         neurons_bv(a).coords=new_tracks{i}.coords;
%         a=a+1;
%     end
% else
%     a=1;
%     neurons_bv=[];
%     for i=1:size(neurons_tmp_unsorted_s,2)
%         neurons_bv(a).coords=short_neurons_tmp_unsorted_s(i).coords;
%         a=a+1;
%     end
% end
% %