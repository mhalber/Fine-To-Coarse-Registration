function [pair, output_matches] = align2view(data, frameID_i, frameID_j, max_dist )
error3D_threshold = 0.05;
error3D_threshold2 = error3D_threshold^2;

%% load two images and depths
image_i = imread(data.image{frameID_i});
image_j = imread(data.image{frameID_j});
scale_x = data.depth_resolution(1)/data.color_resolution(1);
scale_y = data.depth_resolution(2)/data.color_resolution(2);

if ( scale_x ~= 1.0 )
  new_color_res = [data.color_resolution(2) * scale_y, data.color_resolution(1) * scale_x];
	image_i = imresize(image_i, new_color_res);
	image_j = imresize(image_j, new_color_res);
end;

%% convert images to pointclouds. Check if we need to bitshift.
bit_shift = strcmp(data.dataset, 'sun3d');
XYZcam_i = depth2XYZcamera(data.K, depthRead(data.depth{frameID_i}, bit_shift), data.depth_resolution, max_dist);
XYZcam_j = depth2XYZcamera(data.K, depthRead(data.depth{frameID_j}, bit_shift), data.depth_resolution, max_dist);

%% compute SIFT keypoints
[SIFTloc_i,SIFTdes_i] = vl_sift(single(rgb2gray(image_i))) ;
SIFTloc_i = SIFTloc_i([2,1],:);
[SIFTloc_j,SIFTdes_j] = vl_sift(single(rgb2gray(image_j))) ;
SIFTloc_j = SIFTloc_j([2,1],:);

%% SIFT matching
[matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j);

% minNeighboringFrame = 50; % used to be 3, hack to bypass matlab error
minNeighboringMatching = 50;

  if length(matchPointsID_i)<minNeighboringMatching
      fprintf('\nframe %d + %d: too few matching (%d) => relax SIFT threhsold to 0.55\n ', frameID_i, frameID_j , length(matchPointsID_i));
      [matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j, 0.55^2);
      if length(matchPointsID_i)<minNeighboringMatching
        fprintf('with %d matching => relax SIFT threhsold to 0.6\n ', length(matchPointsID_i));
        [matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j, 0.6^2);
        if length(matchPointsID_i)<minNeighboringMatching
          fprintf('with %d matching => relax SIFT threhsold to 0.65\n ', length(matchPointsID_i));
          [matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j, 0.65^2);
          if length(matchPointsID_i)<minNeighboringMatching
            fprintf('with %d matching => relax SIFT threhsold to 0.7\n ', length(matchPointsID_i));
            [matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j, 0.7^2);
            if length(matchPointsID_i)<minNeighboringMatching
              fprintf('with %d matching => relax SIFT threhsold to 0.75\n ', length(matchPointsID_i));
              [matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j, 0.75^2);
              if length(matchPointsID_i)<minNeighboringMatching
                fprintf('with %d matching => relax SIFT threhsold to 0.8\n ', length(matchPointsID_i));
                [matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(SIFTdes_i, SIFTdes_j, 0.8^2);
              end
            end
          end
        end
      end
      fprintf('with %d matching \n', length(matchPointsID_i));
  end

SIFTloc_i = SIFTloc_i(:,matchPointsID_i);
SIFTloc_j = SIFTloc_j(:,matchPointsID_j);

posSIFT_i = round(SIFTloc_i);
valid_i = (1<=posSIFT_i(1,:)) & (posSIFT_i(1,:)<=size(image_i,1)) & (1<=posSIFT_i(2,:)) & (posSIFT_i(2,:)<=size(image_i,2));
posSIFT_j = round(SIFTloc_j);
valid_j = (1<=posSIFT_j(1,:)) & (posSIFT_j(1,:)<=size(image_i,1)) & (1<=posSIFT_j(2,:)) & (posSIFT_j(2,:)<=size(image_i,2));
valid = valid_i & valid_j;

posSIFT_i = posSIFT_i(:,valid);
SIFTloc_i = SIFTloc_i(:,valid);
posSIFT_j = posSIFT_j(:,valid);
SIFTloc_j = SIFTloc_j(:,valid);

Xcam_i = XYZcam_i(:,:,1);
Ycam_i = XYZcam_i(:,:,2);
Zcam_i = XYZcam_i(:,:,3);
validM_i = logical(XYZcam_i(:,:,4));
ind_i = sub2ind([size(image_i,1) size(image_i,2)],posSIFT_i(1,:),posSIFT_i(2,:));
valid_i = validM_i(ind_i);

Xcam_j = XYZcam_j(:,:,1);
Ycam_j = XYZcam_j(:,:,2);
Zcam_j = XYZcam_j(:,:,3);
validM_j = logical(XYZcam_j(:,:,4));
ind_j = sub2ind([size(image_i,1) size(image_i,2)],posSIFT_j(1,:),posSIFT_j(2,:));
valid_j = validM_j(ind_j);
valid = valid_i & valid_j;

ind_i = ind_i(valid);
P3D_i = [Xcam_i(ind_i); Ycam_i(ind_i); Zcam_i(ind_i)];
ind_j = ind_j(valid);
P3D_j = [Xcam_j(ind_j); Ycam_j(ind_j); Zcam_j(ind_j)];

SIFTloc_i = SIFTloc_i(:,valid);
SIFTloc_j = SIFTloc_j(:,valid);
fprintf('Valid SIFT loc %d - %d', size(SIFTloc_i,2), size(SIFTloc_j,2));
%% align RANSAC
output_matches=[-1, -1, -1, -1];
poor_quality = 0;
try
    [RtRANSAC, inliers] = ransacfitRt([P3D_i; P3D_j], error3D_threshold, 0);
    fprintf('frame %6d + %6d: # ransac inliers = %d/%d = %f%%\n', frameID_i, frameID_j, length(inliers), size(P3D_i,2), length(inliers)/size(P3D_i,2)*100);
    output_matches = [ frameID_i, frameID_j, length(inliers), size(P3D_i,2) ];
    if length(inliers) <= 20
      poor_quality = 1;
    end;
catch
    fprintf('frame %6d + %6d: # ransac FAILURE = 0/%d = %f%%\n\n', frameID_i, frameID_j, size(P3D_i,2), 0);
   
    RtRANSAC = [eye(3) zeros(3,1)];
    pair.Rt = RtRANSAC;
    pair.matches = [SIFTloc_i([2 1],:);P3D_i;SIFTloc_j([2 1],:);P3D_j];
    pair.i = frameID_i;
    pair.j = frameID_j;
    pair.valid = 0;
    fprintf('                      Pair %6d + %6d is invalid!\n', frameID_i, frameID_j );
    return;
end

if ~isempty(P3D_i) && ~isempty(P3D_j)
    valid = (sum((P3D_i - transformRT(P3D_j, RtRANSAC, false)).^2,1) < error3D_threshold2);     % Indices of inlying points

%     figure,
%     imshow(image_i);
%     hold on
%     plot(SIFTloc_i(2,valid),SIFTloc_i(1,valid),'+g', 'MarkerSize',10);
%     plot(SIFTloc_i(2,~valid),SIFTloc_i(1,~valid),'+r');
%     drawnow
%     figure,
%     imshow(image_j);
%     hold on
%     plot(SIFTloc_j(2,valid),SIFTloc_j(1,valid),'+g', 'MarkerSize',10);
%     plot(SIFTloc_j(2,~valid),SIFTloc_j(1,~valid),'+r');
%     drawnow

    size_init = size(P3D_i,2);
    validity = 0;
    if size(P3D_i,2) >= 10
      validity = 1;
    end
    
    if poor_quality
      validity = 0;
    end
    
    SIFTloc_i = SIFTloc_i(:,valid);
    P3D_i     = P3D_i    (:,valid);
    SIFTloc_j = SIFTloc_j(:,valid);
    P3D_j     = P3D_j    (:,valid);
    fprintf('                      # point matches  = %d/%d = %f%%\n', size(P3D_i,2), size_init, size(P3D_i,2) / size_init*100);
    fprintf('\t\t %10.6f %10.6f %10.6f %10.6f\n\t\t %10.6f %10.6f %10.6f %10.6f \n\t\t %10.6f %10.6f %10.6f %10.6f\n\n', RtRANSAC' );

    pair.Rt = RtRANSAC;
    pair.matches = [SIFTloc_i([2 1],:);P3D_i;SIFTloc_j([2 1],:);P3D_j];
    pair.i = frameID_i;
    pair.j = frameID_j;
    pair.valid = validity;
    if ~validity 
      fprintf('                      Pair %6d + %6d is invalid!(Too few inliers)\n', frameID_i, frameID_j );
    end

end



