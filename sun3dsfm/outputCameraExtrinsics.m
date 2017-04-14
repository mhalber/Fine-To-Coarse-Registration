function outputCameraExtrinsics(base_path, name, cameraRtC2W)

cameraRtC2W = permute(cameraRtC2W,[2 1 3]);
cameraRtC2W= reshape(cameraRtC2W,4,[]);
cameraRtC2W = cameraRtC2W';
outTo = fullfile(base_path, 'extrinsics', name);
if ~exist(outTo,'dir')
    mkdir(outTo);
end
fp = fopen(fullfile(outTo,[name '.txt']),'w');
for rowID=1:size(cameraRtC2W,1)
    fprintf(fp, '%.100g %.100g %.100g %.100g\n',cameraRtC2W(rowID,:));
end
fclose(fp);    
disp(fullfile(outTo,[name '.txt']));