function [name_vs_id,directions, guess, ruler]=automatic_identity(neurons,chosen_point,voxel_size,fixed_time,final_time)
name_vs_id = [];

% variables initialization
y_scale_factor = voxel_size(2)/voxel_size(1);
z_scale_factor = voxel_size(3)/voxel_size(1);
chosen_point=chosen_point([2 1]);

% signal correlation matrix and coordinates matrix creation
C = correlation_matrix_gf(neurons,final_time);
nen = size(C,1);
coords = zeros(nen,3);
check = 0;
while check == 0
    check = 1;
    for i=1:nen
        coords(i,:) = [neurons(i).coords(fixed_time,1) neurons(i).coords(fixed_time,2).*y_scale_factor neurons(i).coords(fixed_time,3).*z_scale_factor];
        if sum(isnan(coords(i,:))) > 0
            check = 0;
            fixed_time = fixed_time + 1;
        end
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
[guess(1,:), combinations]=find_bpn1(neurons,coords,C,final_time);
guess(2,:)=find_bpn2_8(neurons,coords,C,guess(1,:),combinations,final_time);

% visual check ------------------------------------------------------------
figure
[planeCoeff, normal_vector] = fitPlaneAndPlot_2(coords(guess(:),:));
hold on
% visualCheck(coords,midpoint,guess,coeff,final_time);
for i = 1:size(coords,1)
    posizione(i) = verificaPosizionePunto_2(coords(i,:)', planeCoeff);
    hold on
    if posizione(i)>0
        plot3(coords(i,1),coords(i,2),coords(i,3),'.b','MarkerSize',12)
    else
        plot3(coords(i,1),coords(i,2),coords(i,3),'.g','MarkerSize',12)
    end
    hold on
end
axis equal

% plot_geometry -----------------------------------------------------------

% projection
midpoint=mean(coords(guess(:),:),1)';
ruler=mean(pdist(coords(guess(:),:)));

% insert axes
clear pr_coords
[coeff, ~, ~] = pca(coords);
[anterdir,posterdir,dorsaldir,ventraldir,leftdir,rightdir]=apdvlr_axes_2(coords,guess,100,1, normal_vector, chosen_point);
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
        cost_matrix(j,i)=pdist([pr_coords(fwpn_candidates(i),:).*[1 2 1]; pr_atlas_coords(fpn_ids(j),:)]);
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
a=size(f_pairing,1);
for i=1:size(b_pairing,1)
    name_vs_id{a,1}=bwpn_candidates(b_pairing(i,2));
    name_vs_id{a,2}={namesstr{bpn_ids(b_pairing(i,1))}};
    a=a+1;
end

directions=[anterdir, dorsaldir, leftdir];


function visualCheck(coords,midpoint,guess,coeff,final_time)
figure;
subplot(2,1,1)
plot3(coords(:,1), coords(:,2), coords(:,3), '.','MarkerSize',24);
hold on;
axis equal

guessed=guess(2,:);
guess=guess(1,:);

% Plot the PC1 axis
coeff = coeff;
pc1_vector = coeff(1,:); % First principal component
pc1_start = midpoint - 150 * pc1_vector'; % Starting point for PC1 line
pc1_end = midpoint + 150 * pc1_vector'; % Ending point for PC1 line
plot3([pc1_start(1), pc1_end(1)], [pc1_start(2), pc1_end(2)], [pc1_start(3), pc1_end(3)], 'r-', 'LineWidth', 2);

% Plot midpoint
plot3(midpoint(1),midpoint(2),midpoint(3),'.','MarkerSize',24,'Color',[0.8 0.7 0.1])

% Plot best guess
best_combo = coords(guess,:);
plot3(coords(guess,1), coords(guess,2), coords(guess,3),'.','MarkerSize',24,'Color','r');

% Plot previous guess
best_combo = coords(guessed,:);
for i_bc = 1:size(best_combo,1)
    dist_proj = projectOnPCA(coeff,best_combo(i_bc,:)',midpoint);
    dist_point_3D = midpoint + dist_proj * pc1_vector';
    plot3(dist_point_3D(1),dist_point_3D(2),dist_point_3D(3),'.','MarkerSize',24,'Color',[0.5 0.9 0.8])
end
plot3(coords(guessed,1), coords(guessed,2), coords(guessed,3),'.','MarkerSize',24,'Color','g');
plot3(mean(coords(:,1)), mean(coords(:,2)), mean(coords(:,3)),'.','MarkerSize',24,'Color','m');

hold off; % Release the hold
xlabel('X'); ylabel('Y'); zlabel('Z'); % Label axes
title('3D Plot with PC1 Axis and Projected Points'); % Add title

subplot(2,1,2)
for i_c=guess
   plot(zscore(sv_var(i_c).coords(1:final_time,5)));
   hold on
end
for i_c=guessed
   plot(zscore(sv_var(i_c).coords(1:final_time,5)));
   hold on
end
ylabel(num2str(guess))
midpoint=mean(coords([guess(:) guessed(:)],:),1)'; % per la rappresentazione devo levarlo
subplot(211);
hold on
plot3(midpoint(1),midpoint(2),midpoint(3),'.y','MarkerSize',15)

