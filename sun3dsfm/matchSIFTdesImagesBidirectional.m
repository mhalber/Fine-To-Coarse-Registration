function [matchPointsID_i, matchPointsID_j] = matchSIFTdesImagesBidirectional(X_i, X_j, distRatio)

if ~exist('distRatio','var')
    distRatio = 0.5^2;
end

% l2 normalize sift features
X_i = normc(single(X_i)); 
X_j = normc(single(X_j)); 

numi = size(X_i,2);
numj = size(X_j,2);

% Modification to RootSIFT: https://www.robots.ox.ac.uk/~vgg/publications/2012/Arandjelovic12/arandjelovic12.pdf
for i = 1:numi
   X_i(:,i) = X_i(:,i) ./ sum(abs(X_i(:,i)));
   X_i(:,i) = sqrt( X_i(:,i) );
end

for j = 1:numj
   X_j(:,j) = X_j(:,j) ./ sum(abs(X_j(:,j)));
   X_j(:,j) = sqrt( X_j(:,j) );
end


kdtree_j = vl_kdtreebuild(X_j);
kdtree_i = vl_kdtreebuild(X_i);

matchPointsID_j = zeros(1,numi);
for i = 1 : numi
    [min_idx_j, min_val_j] = vl_kdtreequery(kdtree_j, X_j, X_i(:,i), 'NumNeighbors', 2);
    if (min_val_j(1) < distRatio * min_val_j(2))
        [min_idx_i, min_val_i] = vl_kdtreequery(kdtree_i, X_i, X_j(:,min_idx_j(1)), 'NumNeighbors', 2);
        if min_idx_i(1) == i &&  min_val_i(1) < distRatio * min_val_i(2)
            matchPointsID_j(i) = min_idx_j(1);
        end
    end
end

valid = (matchPointsID_j~=0);

pointsID_i = 1:numi;
matchPointsID_i = pointsID_i(valid);
matchPointsID_j = matchPointsID_j(valid);


