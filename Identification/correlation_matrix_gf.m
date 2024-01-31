function C=correlation_matrix_gf(n,final_time)
nen=numel(n);
signals{nen,nen}=0;
for i=1:nen
    signal = n(i).coords(100:final_time,5);
    signal = signal - imgaussfilt(signal,80);
    signals{i} = signal;
end

for i=1:nen
    for j=i+1:nen
        signal1 = signals{i};
        signal2 = signals{j};
        valid_times=intersect(n(i).coords(:,4),n(j).coords(:,4));
        crrltn=corrcoef(signal1,signal2);
        C(i,j)=crrltn(2);
        C(j,i)=crrltn(2);
    end
end
