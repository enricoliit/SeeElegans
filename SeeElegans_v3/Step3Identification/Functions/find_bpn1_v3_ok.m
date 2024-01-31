function guess=find_bpn1(sv_var,coords,C)
% std of the whole trace
for i=1:size(coords,1)
    stds(i)=std(sv_var(i).coords(:,5));
end
[~, sorted_ids]=sort(stds,'descend');

% sort according to highest stds
guess=[];
selection=sorted_ids(1:60);
num_neu=3;
additional_factor=1;

% consider all possible combinations
selection(ismember(guess(:),selection))=[];
combinations=nchoosek(selection,num_neu);
qscore=zeros(1,size(combinations,1));
[coeff, ~, ~] = pca(coords);
PLs = zeros(size(combinations));
for i=1:size(combinations,1)
    qdr=combinations(i,:);
    
    % filter out all triplets where either one pair has negative correlation
    if sum(sum(C(qdr,qdr)<0))>0 
        continue
    end
    
    % check distances among triplet elements
    xyz=coords(qdr,:)';
    distances=pdist(xyz');
    sorted_distances=sort(distances)-[30 40 60];
    
    % check projection on principal axis
    projLengths = zeros(1,3);
    a = 1;
    for i_qdr=1:numel(qdr)
        for j_qdr=i_qdr+1:numel(qdr)
            projLengths(a) = projectOnPCA(coeff,coords(qdr(i_qdr),:),coords(qdr(j_qdr),:));
            a = a + 1;
        end
    end
    projLengths = sort(abs(projLengths));
    PLs(i,:) = projLengths;
    distcheck2 = sum(abs(projLengths - [30 35 65]));
    
    %distcheck=std(sorted_distances(1:3));
    distcheck=mean(abs(sorted_distances(1:3)));
    Cs=C(qdr,qdr); Cs(Cs==0)=[];
    
    % calculate score of the triplet AVA, AVE, RIM
    std_multiplier=1;
    for k=combinations(i,:)
        std_multiplier=std_multiplier+std(sv_var(k).coords(:,5));
    end
    qscore(i)=std_multiplier*prod(Cs)./(collinearity_check(xyz)*distcheck2*distcheck*additional_factor);
end
[i1,i2]=max(qscore);
guess=combinations(i2,:);
disp(num2str(i1))
disp(num2str(i2))
