function XYZcamera = depth2XYZcamera(K, depth, resolution, max_dist) 
    [x,y] = meshgrid(1:resolution(1), 1:resolution(2));
    XYZcamera(:,:,1) = (x-K(1,3)).*depth/K(1,1);
    XYZcamera(:,:,2) = (y-K(2,3)).*depth/K(2,2);
    XYZcamera(:,:,3) = depth;
    XYZcamera(:,:,4) = depth~=0 & depth < max_dist;
end
