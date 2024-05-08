function scatter_plot=plotspots(ax, scatter_plot, spots, sizes, z, zratio, factor, colrs)
if nargin < 8
colrs=[1 0 0];
end
if ~isempty(spots)
R = sizes(1).^2 - ((z-spots(:,3)).*zratio).^2;
R(R <= 0) = 0.01;
R = sqrt(R);
hold(ax,'on');
scatter_plot=scatter(ax,spots(:,2),spots(:,1),R.*factor,colrs,'LineWidth',factor/20);
hold(ax,'off');
end