function guess=find_bpn2_7(sv_var,coords,C,guessed,combinations,final_time)
top_num = 20;
num_neu = 3;
selection = zeros(1,numel(guessed).*top_num + 1);
for i = guessed
   [~, corrRanking] = sort(C(i,:),'descend'); 
   
   selection = [selection, corrRanking(1:top_num)];
end
selection = unique(selection);
selection(selection==0)=[];
selection(ismember(selection,guessed))=[];

combinations = nchoosek(selection,num_neu);



% initialize variables
additional_factor=1;
qscore=zeros(1,size(combinations,1));
scores=zeros(size(combinations,1),6);
PLs = zeros(size(combinations));

% Calculate midpoint and relative positions of guessed spots along PC1
midpoint = mean(coords,1);
[coeff, ~, ~] = pca(coords);

% Check distances between first guessed spots
xyz_g=coords(guessed,:)';
distances_guessed=abs(pdist(xyz_g'));
sorted_distances_guessed=sort(distances_guessed);

% Check projections on PC1 of distances between first guessed spots
projLengths_guessed = zeros(1,3);
a = 1;
for i_qdr=1:numel(guessed)
    for j_qdr=i_qdr+1:numel(guessed)
        projLengths_guessed(a) = projectOnPCA(coeff,coords(guessed(i_qdr),:),coords(guessed(j_qdr),:));
        a = a + 1;
    end
end
projLengths_guessed = (abs(projLengths_guessed));
sorted_projLengths_guessed = sort(projLengths_guessed);

% Check projections on PC1 of first guessed spots
coords_guessed = coords(guessed,:);
dist_proj_guessed = zeros(1,numel(guessed));
for i = 1:size(coords_guessed,1)
    dist_proj_guessed(i) = projectOnPCA(coeff,coords_guessed(i,:),midpoint);
end
sorted_dist_proj_guessed=sort(dist_proj_guessed);

Cg = sum(sum(C(guessed,guessed)));

% Loop over all combinations
for i=1:size(combinations,1)
    qdr=combinations(i,:);

    
    % Filter out all triplets where either one pair has negative correlation
    if sum(sum(C(qdr,qdr)<0))>0 
        continue
    end
    if sum(sum(C(qdr,qdr))) < 0.8 * Cg
        continue
    end

    
    
    
    % Check distances among triplet elements
    xyz=coords(qdr,:)';
    distances=abs(pdist(xyz'));
    sorted_distances = sort(distances);
    if any(sorted_distances < 5)
       continue 
    end
    sorted_distances_check = sorted_distances - sorted_distances_guessed;
    
    
    
    
    % Check distances on projection on principal axis
    projLengths = zeros(1,3);
    a = 1;
    for i_qdr=1:numel(qdr)
        for j_qdr=i_qdr+1:numel(qdr)
            projLengths(a) = projectOnPCA(coeff,coords(qdr(i_qdr),:),coords(qdr(j_qdr),:));
            a = a + 1;
        end
    end
    projLengths = (abs(projLengths));
    sorted_projLengths = sort(projLengths);
    PLs(i,:) = projLengths;
    
    
    
    % Filter out any triplet that's not arranged along the PC1 of the spots
    [~, im2] = max(distances);
    if projLengths(im2)/distances(im2) < 0.85
       continue 
    end
    distcheck2 = sum(abs(sorted_projLengths - sorted_projLengths_guessed));
    
    
    
    
    % Calculate relative positions along PC1 for current triplet
    coords_qdr = coords(qdr,:);
    dist_proj_qdr = zeros(1,numel(guessed));
    for i_qdr = 1:numel(qdr)
        dist_proj_qdr(i_qdr) = projectOnPCA(coeff,coords_qdr(i_qdr,:),midpoint);
    end
    sorted_dist_proj_qdr=sort(dist_proj_qdr);
    distchecks3 = (sorted_dist_proj_guessed - sorted_dist_proj_qdr).^2;
    RLs(i,:) = sorted_dist_proj_qdr;

%     if any(distchecks3>1200)
%         continue
%     end
    
    distcheck3 = sum(distchecks3);
    
    distcheck=mean(abs(sorted_distances_check(1:3)));
    Cs=C(qdr,qdr); Cs(Cs==0)=[];
    
    % calculate score of the triplet AVA, AVE, RIM/AIB
    std_multiplier=1;
    for k=combinations(i,:)
        std_multiplier=std_multiplier+std(sv_var(k).coords(:,5),'omitnan');
    end
    qscore(i)=std_multiplier*prod(Cs)./((collinearity_check(xyz)*distcheck3*distcheck2*distcheck*additional_factor));
    scores(i,:) = [std_multiplier prod(Cs) collinearity_check(xyz) distcheck3 distcheck2 distcheck];
end
[i1,i2]=max(qscore);
[i1s,i2s]=sort(qscore); 
% combinations(i2s(~isnan(i1s)),:)
guess=combinations(i2,:);

% Plotting the original points
figure;
subplot(2,1,1)
plot3(coords(:,1), coords(:,2), coords(:,3), '.','MarkerSize',24);
hold on;
axis equal
% Plot the PC1 axis
pc1_vector = coeff(:,1); % First principal component
pc1_start = midpoint - 150 * pc1_vector'; % Starting point for PC1 line
pc1_end = midpoint + 150 * pc1_vector'; % Ending point for PC1 line
plot3([pc1_start(1), pc1_end(1)], [pc1_start(2), pc1_end(2)], [pc1_start(3), pc1_end(3)], 'r-', 'LineWidth', 2);
% Plot midpoint
plot3(midpoint(1),midpoint(2),midpoint(3),'.','MarkerSize',24,'Color',[0.8 0.7 0.1])
% Plot best guess
best_combo = coords(guess,:);
for i = 1:size(best_combo,1)
    dist_proj = projectOnPCA(coeff,best_combo(i,:),midpoint);
    dist_point_3D = midpoint + dist_proj * pc1_vector';
    plot3(dist_point_3D(1),dist_point_3D(2),dist_point_3D(3),'.','MarkerSize',24,'Color',[0.9 0.5 0.8])
end
plot3(coords(guess,1), coords(guess,2), coords(guess,3),'.','MarkerSize',24,'Color','r');
% Plot previous guess
best_combo = coords(guessed,:);
for i = 1:size(best_combo,1)
    dist_proj = projectOnPCA(coeff,best_combo(i,:),midpoint);
    dist_point_3D = midpoint + dist_proj * pc1_vector';
    plot3(dist_point_3D(1),dist_point_3D(2),dist_point_3D(3),'.','MarkerSize',24,'Color',[0.5 0.9 0.8])
end
plot3(coords(guessed,1), coords(guessed,2), coords(guessed,3),'.','MarkerSize',24,'Color','g');
plot3(mean(coords(:,1)), mean(coords(:,2)), mean(coords(:,3)),'.','MarkerSize',24,'Color','m');

hold off; % Release the hold
xlabel('X'); ylabel('Y'); zlabel('Z'); % Label axes
title('3D Plot with PC1 Axis and Projected Points'); % Add title

subplot(2,1,2)
for i=guess
   plot(zscore(sv_var(i).coords(1:final_time,5)));
   hold on
end
for i=guessed
   plot(zscore(sv_var(i).coords(1:final_time,5)));
   hold on
end
ylabel(num2str(guess))
midpoint=mean(coords([guess(:) guessed(:)],:),1)'; % per la rappresentazione devo levarlo
subplot(211);
hold on
plot3(midpoint(1),midpoint(2),midpoint(3),'.y','MarkerSize',15)
