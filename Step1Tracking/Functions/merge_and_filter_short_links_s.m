function filtered_exptracelist=merge_and_filter_short_links_s(tracklist,segmentlists_time,segmentlists,length_filter,tmax_local)

fet=1;
filtered_exptracelist=[];
for tracklist_id=1:numel(tracklist)
    new_track=NaN.*ones(1,tmax_local);
    for segment_id=tracklist{tracklist_id}
        start_time=(segmentlists_time{segment_id}(1));
        end_time=(segmentlists_time{segment_id}(end));
        new_track(start_time:end_time)=segmentlists{segment_id};
    end
    if sum(~isnan(new_track))<length_filter
        continue
    end    
    filtered_exptracelist{fet}=new_track;
    fet=fet+1;
end
