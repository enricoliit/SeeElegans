function C=xcorrelation_matrix(n)
nen=numel(n);
for i=1:nen
    for j=i+1:nen
        valid_times=intersect(n(i).coords(:,4),n(j).coords(:,4));
        crrltn_p=max(xcorr(n(i).coords(valid_times,5),n(j).coords(valid_times,5)));
        crrltn_n=max(xcorr(n(i).coords(valid_times,5),-1.*n(j).coords(valid_times,5)));
        [crrSign, signPos] = max([crrltn_p crrltn_n]);
        if signPos == 1
            C(i,j)=crrltn_p;
            C(j,i)=crrltn_p;
        else
            C(i,j)=-1*crrltn_n;
            C(j,i)=-1*crrltn_n;
        end
    end
end
