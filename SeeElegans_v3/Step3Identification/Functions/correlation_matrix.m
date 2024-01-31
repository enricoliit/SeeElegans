function C=correlation_matrix(n)
nen=numel(n);
for i=1:nen
    for j=i+1:nen
        valid_times=intersect(n(i).coords(:,4),n(j).coords(:,4));
        crrltn=corrcoef(n(i).coords(valid_times,5),n(j).coords(valid_times,5));
        C(i,j)=crrltn(2);
        C(j,i)=crrltn(2);
    end
end
