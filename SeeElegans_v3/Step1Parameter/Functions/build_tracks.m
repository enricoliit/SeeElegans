function tracklist=build_tracks(segment_links)

origin_id=1;
tracklist{origin_id}=[];

while ~isempty(segment_links)
    current_index=segment_links(1,1);
    target_index=segment_links(1,2);
    tracklist{origin_id}=[current_index, target_index];
    segment_links(1,:)=[];
    next_segment=find(segment_links(:,1)==target_index);
    while ~isempty(next_segment)
        current_index=segment_links(next_segment,1);
        target_index=segment_links(next_segment,2);
        tracklist{origin_id}=[tracklist{origin_id}, target_index];
        segment_links(next_segment,:)=[];
        next_segment=find(segment_links(:,1)==target_index);
    end 
    tracklist{origin_id}(tracklist{origin_id}==0)=[];
    origin_id=origin_id+1;
end