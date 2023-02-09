function guess=find_bpn1(sv_var,coords,C)
for i=1:size(coords,1)
    stds(i)=std(sv_var(i).coords(:,5));
end
[~, sorted_ids]=sort(stds,'descend');

guess=[];
selection=sorted_ids(1:50);
num_neu=3;
additional_factor=1;

selection(ismember(guess(:),selection))=[];
combinations=nchoosek(selection,num_neu);
qscore=zeros(1,size(combinations,1));
for i=1:size(combinations,1)
    qdr=combinations(i,:);
    if sum(sum(C(qdr,qdr)<0))>0
        continue
    end
    xyz=coords(qdr,:)';
    distances=pdist(xyz');
    sorted_distances=sort(distances);
    distcheck=std(sorted_distances(1:3));
    Cs=C(qdr,qdr); Cs(Cs==0)=[];
    std_multiplier=1;
    for k=combinations(i,:)
        std_multiplier=std_multiplier*std(sv_var(k).coords(:,5));
    end
    qscore(i)=std_multiplier*prod(Cs)./(collinearity_check(xyz)*distcheck*additional_factor);
end
[~,i2]=max(qscore);
guess=combinations(i2,:);
