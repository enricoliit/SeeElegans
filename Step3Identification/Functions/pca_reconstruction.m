function varb1=pca_reconstruction(varb,rank)
if nargin<2
   rank=4; 
end

time1=max(varb(1).coords(:,4));
for i=2:numel(varb)
    time1=min([time1 max(varb(i).coords(:,4))]);
end
raster=zeros(numel(varb),time1,3);
reraster=zeros(numel(varb),time1,3);
for i=1:numel(varb)
    raster(i,:,1)=varb(i).coords(1:time1,1);
    raster(i,:,2)=varb(i).coords(1:time1,2);
    raster(i,:,3)=varb(i).coords(1:time1,3);
end
for i=1:3
factors = raster(:,:,i);
observations = factors ;
[coeff, score, ~, ~, ~, mu] = pca(observations);
reconstructed = score(:,1:rank) * coeff(:,1:rank)' + repmat(mu, size(factors,1), 1);
sum((observations - reconstructed).^2)
reraster(:,:,i)=reconstructed;
end
varb1=varb;
for i=1:numel(varb)
    varb1(i).coords(1:time1,1)=reraster(i,:,1);
    varb1(i).coords(1:time1,2)=reraster(i,:,2);
    varb1(i).coords(1:time1,3)=reraster(i,:,3);
end