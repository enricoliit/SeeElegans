function spots=bd3d
global framenum
global sizes
global thresholds
global sigmas
global wholestack
global zsizes

A=wholestack(:,:,:,framenum);
h=-fspecial3('log',[sizes sizes zsizes],sigmas);
If=imfilter(A,h,'symmetric');
regional_max_frame=imregionalmax(If,26);
try
    [x_i,y_i,z_i]=ind2sub(size(A),find(regional_max_frame.*A>thresholds));
    spots_unfiltered=[x_i,y_i,z_i];
    max_in_A_times_If=regional_max_frame.*A;
    linear_maxima=find(max_in_A_times_If>thresholds);
    q_i=If(linear_maxima);
    qI_i=q_i.*A(linear_maxima);
catch
    nat=whos('A');
    class_of_frame=nat.class;
    eval(['regional_max_frame=' class_of_frame '(regional_max_frame);']);
    [x_i,y_i,z_i]=ind2sub(size(A),find(regional_max_frame.*A>thresholds));
    spots_unfiltered=[x_i,y_i,z_i];
    max_in_A_times_If=regional_max_frame.*A;
    linear_maxima=find(max_in_A_times_If>thresholds);
    q_i=If(linear_maxima);
    qI_i=q_i.*A(linear_maxima);
end

%% proximity check
overlap_of_neurons=NaN*ones(numel(spots_unfiltered));
a=1;
for i=1:size(spots_unfiltered,1)
   for j=i+1:size(spots_unfiltered,1)
       overlap_of_neurons(i,j)=sqrt( (spots_unfiltered(i,1)-spots_unfiltered(j,1)).^2 + (spots_unfiltered(i,2)-spots_unfiltered(j,2)).^2 );
       if overlap_of_neurons(i,j)<round(sizes./2) && abs(spots_unfiltered(i,3)-spots_unfiltered(j,3))<round(zsizes./2)
           too_close_pairs(a,:)=[i, j];
           a=a+1;
       end
   end
end

%% remove redundant tracks
try
    a=1;
    for i=1:size(too_close_pairs,1)
        [~, minpos]=min([qI_i(too_close_pairs(i,1)) qI_i(too_close_pairs(i,2))]);
        traces_to_be_deleted(a)=too_close_pairs(i,minpos);
        a=a+1;
    end
end

a=1;
if ~exist('traces_to_be_deleted','var')
    traces_to_be_deleted=[];
    spots=spots_unfiltered;
else
    for i=1:size(spots_unfiltered,1)
        if ~ismember(i,traces_to_be_deleted)
            spots(a,:)=spots_unfiltered(i,:);
            a=a+1;
        end
    end
end















