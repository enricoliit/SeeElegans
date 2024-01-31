function [D, dd, xyz00]=collinearity_check(xyz)
xyz0 = mean(xyz,2);
A = xyz-xyz0;
[U,~,~] = svd(A);
d = U(:,1);
t = d'*A;
t1 = min(t);
t2 = max(t);
xzyl = xyz0 + [t1,t2].*d;
X1=xzyl(:,1);
X2=xzyl(:,2);
D=0;
for i=1:size(xyz,2)
    X0=xyz(:,i);
    D=D+sqrt(sum(cross((X0-X1),(X0-X2)).^2)/sum((X2-X1).^2));
end
if nargout>1
    dd=d;
    xyz00=xyz0;
end