function guess=find_bpn2_2(neurons,coords,C,guessed)
for i=1:size(coords,1)
    stds(i)=std(neurons(i).coords(:,5));
end
[~, sorted_ids]=sort(stds,'descend');

guess=[];
selection=sorted_ids(1:50);
num_neu=3;
additional_factor=1;

selection(ismember(guess(:),selection))=[];
combinations=nchoosek(selection,num_neu);
valid_times=neurons(i).coords(1,4);
for i=1:numel(neurons)
    valid_times=intersect(valid_times,neurons(i).coords(:,4));
end
for i=1:numel(neurons)
    tr_signal=neurons(i).coords(valid_times,5);
    tr_signal(isnan(tr_signal))=[];
    raster(i,:)=zscore(tr_signal);
end
kmns=clusterdata(raster,'Distance','correlation','maxclust',4);
S=mdwtcluster(raster,'maxclust',4,'pdist','correlation');
target_cluster=S.IdxCLU(guessed(1,1));
xyz=coords(guessed,:)';
[~,mxyz,dxyz]=collinearity_check(xyz);
slope_factor=0;
[clcxyz,mclcxyz,~]=collinearity_check(xyz);
c_to_remove=[];
for i=1:size(combinations,1)
    qdr=combinations(i,:);
    for j=1:numel(qdr)
        if C(qdr(j),qdr(j))<0
            c_to_remove=[c_to_remove, i];
            break
        end
        if ismember(qdr(j),guessed(1,:))
            c_to_remove=[c_to_remove, i];
            break
        end
        if S.IdxCLU(qdr(j),1)~=target_cluster
            c_to_remove=[c_to_remove, i];
            break
        end
    end
end

combinations(c_to_remove,:)=[];
a=1;
clear qscore qscores
for i=1:size(combinations,1)
    qdr=combinations(i,:);
    xyz=coords(qdr,:)';
    distances=pdist(xyz');
    if sum(distances<10)>0
        continue
    end
    sorted_distances=sort(distances);
    distcheck=std(sorted_distances(1:2));
    Cs=C(qdr,qdr); Cs(Cs==0)=[];
    std_multiplier=1;
    for k=combinations(i,:)
        std_multiplier=std_multiplier*std(neurons(k).coords(:,5));
    end
    [clcxyz,mclcxyz,dclcxyz]=collinearity_check(xyz);
    dist_g1=zeros(1,numel(combinations(i,:)));
    taken=[];
    candidates=guessed(1,:);
    a=1;
    for m=combinations(i,:)
        candidates(ismember(candidates,taken))=[];
        single_dist_g1=(squareform(pdist([coords(candidates,:); coords(m,:)])));
        single_dist_g1(single_dist_g1==0)=NaN;
        [to_remove,single_dist_g1]=min(single_dist_g1(:,end));
        taken=[taken, to_remove];
        dist_g1(a)=single_dist_g1;
        a=a+1;
    end
    additional_factor=std(dist_g1);
    additional_factor_s=std(dist_g1);
    additional_factor_m=mean(dist_g1);
    slope_factor=sum((mxyz(1)-mclcxyz(1)).^2);
    qscore(i)=(sum((dclcxyz-dxyz).^2)+std_multiplier*prod(Cs))./(collinearity_check(xyz)*distcheck*additional_factor*additional_factor_s+slope_factor+additional_factor_m);
    qscores(i,:)=[sum((dclcxyz-dxyz).^2),std_multiplier,prod(Cs),collinearity_check(xyz),distcheck,additional_factor,additional_factor_s,slope_factor,additional_factor_m];
    a=a+1;
end
[~,i2]=max(qscore);
guess=combinations(i2,:);
