function outputConf(conf_name, data, cameraRt2CW)

fid = fopen( conf_name, 'w' );
fprintf(fid, 'dataset %s\n', data.dataset);
fprintf(fid, 'n_images %d\n', data.n_images);
fprintf(fid, 'intrinsics %s\n', data.intrinsics_name);
fprintf(fid, 'color_resolution %d %d\n', data.color_resolution(1), data.color_resolution(2) );
fprintf(fid, 'depth_resolution %d %d\n', data.depth_resolution(1), data.depth_resolution(2) );
fprintf(fid, 'color_directory %s\n', data.image_directory );
fprintf(fid, 'depth_directory %s\n', data.depth_directory );
fprintf(fid, 'depth_directory matches/rgbdsfm_adjacent_matches\n' );
fprintf(fid, '\n' );

for i = 1:data.n_images
    [~,color_name, color_ext] = fileparts(data.image{i});
    [~,depth_name, depth_ext] = fileparts(data.depth{i});
    T = [cameraRt2CW(:,:,i); 0 0 0 1];
    C = [ 1 0 0 0 ; 0 -1 0 0 ; 0 0 -1 0 ; 0 0 0 1];
    T = C * T * C;
    
    fprintf(fid, 'scan %s %s ', strcat(depth_name, depth_ext), strcat(color_name,color_ext) );
    fprintf(fid, '%f %f %f %f ', T(1,:) );
    fprintf(fid, '%f %f %f %f ', T(2,:) );
    fprintf(fid, '%f %f %f %f ', T(3,:) );
    fprintf(fid, '%f %f %f %f\n', T(4,:) );
end

fclose(fid);
fprintf('Wrote conf file to %s\n', conf_name);
end