function links=LAP_linker_local_s(framenum_id,non_linking_cost,spots_s,xy_pixel_size,z_pixel_size)
try
    movement_correction_factor=0.45;
    z_pixel_size=z_pixel_size*movement_correction_factor;

    % Subpixel centroid extraction in previous frame (t-1), (xc_t0,yc_t0,zc_t0)
    xc_t0 = xy_pixel_size .* spots_s{framenum_id-1}(:,1);
    yc_t0 = xy_pixel_size .* spots_s{framenum_id-1}(:,2);
    zc_t0 =  z_pixel_size .* spots_s{framenum_id-1}(:,3);
    
    % Subpixel centroid extraction in previous frame (t1), (xc_t1,yc_t1,zc_t1)
    xc_t1 = xy_pixel_size .* spots_s{framenum_id}(:,1);
    yc_t1 = xy_pixel_size .* spots_s{framenum_id}(:,2);
    zc_t1 =  z_pixel_size .* spots_s{framenum_id}(:,3);
    
    % cost matrix calculation
    C=zeros(length(xc_t1),numel(xc_t0));
    for k = 1:numel(xc_t0)
        for h = 1:length(xc_t1)
            C(h,k) = vecnorm([xc_t0(k), yc_t0(k), zc_t0(k)]-[xc_t1(h), yc_t1(h), zc_t1(h)]);
        end
    end
    C(C > non_linking_cost) = Inf;
    C=C.^2;
    highest_cost=max(C (C ~= Inf));
    UnmatchedCost=double(highest_cost*1.05);

    % lowest cost solution
    UnmatchedCost(isempty(UnmatchedCost))=0;
    links=matchpairs(C,UnmatchedCost);
    links=links(:,[2 1]);

catch

    links=[];

end