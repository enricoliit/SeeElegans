function [anterdir,posterdir,dorsaldir,ventraldir,leftdir,rightdir]=apdvlr_axes(coords,guess,scale_factor,str,ant_dx)
midpoint=mean(coords,1,'omitnan')';
[coefs,~]=pca(coords);
dorsaldir_counter=0;
for i=guess(:)'
    dorsaldir_counter=dorsaldir_counter+dot(midpoint,coefs(:,2))-dot(coefs(:,2),coords(i,:));
end
dorsaldir=1.*sign(dorsaldir_counter).*coefs(:,2);
ventraldir=-1.*sign(dorsaldir_counter).*coefs(:,2);

if ant_dx==1
    anterdir=-1.*sign(dorsaldir_counter).*coefs(:,1);
    posterdir=1.*sign(dorsaldir_counter).*coefs(:,1);
    leftdir=-1.*sign(dorsaldir_counter).*coefs(:,3);
    rightdir=1.*sign(dorsaldir_counter).*coefs(:,3);
else
    anterdir=1.*sign(dorsaldir_counter).*coefs(:,1);
    posterdir=-1.*sign(dorsaldir_counter).*coefs(:,1);
    leftdir=1.*sign(dorsaldir_counter).*coefs(:,3);
    rightdir=-1.*sign(dorsaldir_counter).*coefs(:,3);
end
if str==1
    hold on
    plot3([0, ventraldir(1)]'.*scale_factor',[0, ventraldir(2)].*scale_factor',[0, ventraldir(3)].*scale_factor');
    hold on
    plot3([0, dorsaldir(1)]'.*scale_factor',[0, dorsaldir(2)].*scale_factor',[0, dorsaldir(3)].*scale_factor');
    hold on
    text(ventraldir(1).*scale_factor,ventraldir(2).*scale_factor,ventraldir(3).*scale_factor,'V')
    hold on
    text(dorsaldir(1).*scale_factor,dorsaldir(2).*scale_factor,dorsaldir(3).*scale_factor,'D')
    hold on
    plot3([0, anterdir(1)]'.*scale_factor',[0, anterdir(2)].*scale_factor',[0, anterdir(3)].*scale_factor');
    hold on
    plot3([0, posterdir(1)]'.*scale_factor',[0, posterdir(2)].*scale_factor',[0, posterdir(3)].*scale_factor');
    hold on
    plot3([0, rightdir(1)]'.*scale_factor',[0, rightdir(2)].*scale_factor',[0, rightdir(3)].*scale_factor');
    hold on
    plot3([0, leftdir(1)]'.*scale_factor',[0, leftdir(2)].*scale_factor',[0, leftdir(3)].*scale_factor');
    hold on
    text(anterdir(1).*scale_factor,anterdir(2).*scale_factor,anterdir(3).*scale_factor,'A')
    hold on
    text(posterdir(1).*scale_factor,posterdir(2).*scale_factor,posterdir(3).*scale_factor,'P')
    hold on
    text(leftdir(1).*scale_factor,leftdir(2).*scale_factor,leftdir(3).*scale_factor,'L')
    hold on
    text(rightdir(1).*scale_factor,rightdir(2).*scale_factor,rightdir(3).*scale_factor,'R')
end