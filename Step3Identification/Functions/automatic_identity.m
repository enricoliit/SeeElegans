function [name_vs_id,directions]=automatic_identity(neurons,chosen_point,voxel_size,fixed_time)
name_vs_id=[];

% variables initialization
y_scale_factor=voxel_size(2)/voxel_size(1);
z_scale_factor=voxel_size(3)/voxel_size(1);

% signal correlation matrix and coordinates matrix creation
C=correlation_matrix(neurons);
nen=size(C,1);
coords=zeros(nen,3);
for i=1:nen
    try
        coords(i,:)=[neurons(i).coords(fixed_time,1) neurons(i).coords(fixed_time,2).*y_scale_factor neurons(i).coords(fixed_time,3).*z_scale_factor];
    catch
    end
end

% check head direction
ant_dx=0;
if ~isempty(chosen_point)
    coordsmeanpoint=mean(coords(:,1));
    if chosen_point(1)-coordsmeanpoint>0
        ant_dx=1;
        disp('testa a destra');
    end
end
% atlas initialization
load('./Functions/atlas.mat');
bpn_ids_s=bpn_ids(1:end-2);
midpoint_str=mean(m_coords_str_corrected(bpn_ids_s,:),1);
ruler_str=mean(pdist(m_coords_str_corrected(bpn_ids_s,:)));
for i=1:numel(m_coords_str_corrected(:,1))
    pr_atlas_coords(i,:)=dirprojectcoords(anterdir_str,dorsaldir_str,leftdir_str,midpoint_str,m_coords_str_corrected(i,:))./ruler_str;
end

%% neuron arrangement in recording
% C=correlation_matrix(neurons);
guess(1,:)=find_bpn1(neurons,coords,C);
guess(2,:)=find_bpn2_2(neurons,coords,C,guess(1,:));

% projection
midpoint=mean(coords(guess(:),:),1)';
ruler=mean(pdist(coords(guess(:),:)));
clear pr_coords
[anterdir,posterdir,dorsaldir,ventraldir,leftdir,rightdir]=apdvlr_axes(coords,guess,1,1,ant_dx);
for i=1:numel(coords(:,1))
    pr_coords(i,:)=dirprojectcoords(anterdir,dorsaldir,leftdir,midpoint,coords(i,:))./ruler;
end

%% comparison
identified=guess(:);
id_check=1:numel(neurons);
id_check(ismember(id_check,identified))=[];
for i=id_check
    id_score(i)=sum(C(i,identified));
end
clear cost_matrix
fwpn_candidates=find(id_score<-std(id_score(id_score<0),0,'omitnan'));
cost_matrix=[];
for i=1:numel(fwpn_candidates)
    for j=1:numel(fpn_ids)
        cost_matrix(j,i)=pdist([pr_coords(fwpn_candidates(i),:).*[1 2 1];pr_atlas_coords(fpn_ids(j),:)]);
    end
end
f_pairing=[];
if ~isempty(cost_matrix)
    f_pairing=matchpairs(cost_matrix,2);
end


%%
identified=guess(:);
id_check=1:numel(neurons);
for i=id_check
    id_score(i)=sum(C(i,identified));
end
clear cost_matrix
bwpn_candidates=find(id_score>std(id_score,0,'omitnan'));
cost_matrix=zeros(numel(bpn_ids),numel(bwpn_candidates));
for i=1:numel(bwpn_candidates)
    for j=1:numel(bpn_ids)
        cost_matrix(j,i)=pdist([pr_coords(bwpn_candidates(i),:).*[1 3 1];pr_atlas_coords(bpn_ids(j),:)])./id_score(bwpn_candidates(i));
    end
end

b_pairing=matchpairs(cost_matrix,2);

%% name id pairing
for i=1:size(f_pairing,1)
    name_vs_id{i,1}=fwpn_candidates(f_pairing(i,2));
    name_vs_id{i,2}={namesstr{fpn_ids(f_pairing(i,1))}};
end
a=1;
for i=1:size(b_pairing,1)
    name_vs_id{a,1}=bwpn_candidates(b_pairing(i,2));
    name_vs_id{a,2}={namesstr{bpn_ids(b_pairing(i,1))}};
    a=a+1;
end

directions=[anterdir, dorsaldir, leftdir];