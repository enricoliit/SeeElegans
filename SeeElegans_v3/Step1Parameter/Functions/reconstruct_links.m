function reconstructed_links=reconstruct_links(mutual_position_based_links,max_gap_closing_mutual_distance)

if ~isempty(mutual_position_based_links)
    links_backup=mutual_position_based_links;
    ex_indices = find(mutual_position_based_links(:,3)>max_gap_closing_mutual_distance);
    mutual_position_based_links(ex_indices,:) = [];
    mutual_position_based_links=mutual_position_based_links(any(mutual_position_based_links,2),:);
    clear reconstruction reconstructed_links
    a=0;
    h=waitbar(0,'Rebuilding Tracks...');
    for i=1:size(mutual_position_based_links,1)
        waitbar(i/size(mutual_position_based_links,1),h);
        if mutual_position_based_links(i,1)>0
            a=a+1;
            reconstructed_links{a}=[ mutual_position_based_links(i,1) mutual_position_based_links(i,2)];
            last=mutual_position_based_links(i,2);
            mutual_position_based_links(i,:)=[0 0 0];
            if last==0
                continue
            end
            while last~=0
                found=find(mutual_position_based_links(:,1)==last);
                if numel(found)>0
                    reconstructed_links{a}=[reconstructed_links{a}, mutual_position_based_links(found,2)];
                    last=mutual_position_based_links(found,2);
                    mutual_position_based_links(found,:)=[0 0 0];
                else
                    break
                end
            end
            reconstructed_links{a}=[reconstructed_links{a}, links_backup(found,2)];
        end
    end
    close(h)
end
