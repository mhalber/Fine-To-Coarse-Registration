function outputPly(PLYfilename, data, cameraRtC2W)

tic;
fprintf('Writing ply point cloud file: ');

scale_x = data.depth_resolution(1)/data.color_resolution(1);
scale_y = data.depth_resolution(2)/data.color_resolution(2);

frameInterval = 10;

pointCount = round(min(data.depth_resolution(1)*data.depth_resolution(2)*0.5,1000000/(length(data.image)/frameInterval)));

bitshift = strcmp(data.dataset, 'sun3d');
dataChunk = uint8([]);
for frameID=1:frameInterval:length(data.image)

    image = imread(data.image{frameID});
    if ( scale_x ~= 1.0 || scale_y ~= 1.0 )
        new_color_res = [data.color_resolution(2) * scale_y, data.color_resolution(1) * scale_x];
        image = imresize(image, new_color_res);
    end;
    XYZcam = depth2XYZcamera(data.K, depthRead(data.depth{frameID}, bitshift), data.depth_resolution, 5.0);       
    
    XYZcam = reshape(XYZcam,data.depth_resolution(1)*data.depth_resolution(2),4)';
    isValid = find(XYZcam(4,:));
    XYZcam = transformRT(XYZcam(1:3,:),cameraRtC2W(:,:,frameID));
    try
        selID = randsample(length(isValid),pointCount);
    catch
        selID = randsample(length(isValid),pointCount,true);
    end
    isValid = isValid(selID);
 
    RGB =  reshape(image, data.depth_resolution(1)*data.depth_resolution(2), 3)';
    
    dataChunk = [dataChunk [reshape(typecast(reshape(single(XYZcam(:,isValid)),1,[]),'uint8'),3*4,[]); RGB(:,isValid)]];
end


file = writePLYhead(PLYfilename, size(dataChunk,2));
fwrite(file, dataChunk,'uint8');
fclose(file);

toc;
