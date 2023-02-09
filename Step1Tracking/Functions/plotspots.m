function plotspots(spots)
global ax
global sizes
global I
global text_spot_num
global slice
global zsizes
ax; cla; imagesc(I), axis equal, axis off, hold on

for i=1:zsizes
    selspots=(abs(spots(:,3)-slice)==(i-1));
    visspots=spots(selspots,:);
    viscircles(visspots(:,[2 1]),ones(size(visspots,1),1).*round(sizes)./(2*i));
end
selspots_off_plane=(abs(spots(:,3)-slice)>(zsizes-1));
visspots_off_plane=spots(selspots_off_plane,:);
viscircles(visspots_off_plane(:,[2 1]),ones(size(visspots_off_plane,1),1));

set(text_spot_num,'String',['spots = ' num2str(size(spots,1))]);
end
