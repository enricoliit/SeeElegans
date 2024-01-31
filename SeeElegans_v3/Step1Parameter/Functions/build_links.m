function [segmentlists_time,segmentlists]=build_links(links)

links{end+1}=NaN;
segmentlists=[];
segmentlists_time=[];
origin_id=1;
for i=2:numel(links)                   
    origin_id=1;
    traces_id=numel(segmentlists)+1;
    while ~isempty(links{i}) && sum(~isnan(links{i}(:,1)))~=0
        origin=links{i}(origin_id,1);
        if isnan(origin)
            origin_id=origin_id+1;
            continue
        end
        segmentlists{traces_id}(1)=origin;
        segmentlists_time{traces_id}(1)=i-1;
        target=links{i}(origin_id,2);
        links{i}(origin_id,:)=NaN;
        search_t=i+1;        
        try
            linking=find(links{search_t}(:,1)==target);
        catch
            linking=[];
        end
        while ~isempty(linking) && ~isempty(links{search_t})
            linking=linking(1);
            segmentlists{traces_id}(end+1)=links{search_t}(linking,1);
            segmentlists_time{traces_id}(end+1)=search_t-1;
            origin=links{search_t}(linking,1);
            target=links{search_t}(linking,2);
            links{search_t}(linking,:)=NaN;
            search_t=search_t+1;
            try
                linking=find(links{search_t}(:,1)==target);
            end
        end
        segmentlists{traces_id}(end+1)=target;
        segmentlists_time{traces_id}(end+1)=search_t-1;
        origin_id=origin_id+1;
        traces_id=traces_id+1;
    end
end
