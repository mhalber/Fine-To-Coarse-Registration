function cameraRtC2W=sun3Dsfm( conf_name, base_path, sun3dsfm_path )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conf_name - name of the configuration file you wish to load
% base_path - directory path to a folder where <conf_name> is located
% sun3dsfm_path - path to where this file is located.

sequence_tokens = strsplit( base_path, '/');
n_tokens = size( sequence_tokens, 2 );
sequence_name = sequence_tokens{n_tokens};
if isempty(sequence_name)
    sequence_name = sequence_tokens{n_tokens-1};
end

if  ~exist('conf_name', 'var')
    conf_name = strcat( sequence_name ,'.conf' );
end;

if ~exist('base_path', 'var')
    base_path = pwd;
end;

if  ~exist('sun3dsfm_path', 'var')
    % Note you might wish to enter your path here.
    sun3dsfm_path='';
end;

if ~exist('frameIDs','var')
    frameIDs = [];
end

if ~exist('doLoopClosureMatches','var')
    doLoopClosureMatches = 0;
end

if ~exist('doLongTrack','var')
    doLongTrack = 0;
end

%% library setup	
addpath(genpath(fullfile(sun3dsfm_path,'lib/estimateRigidTransform')));	
addpath(genpath(fullfile(sun3dsfm_path,'lib/peter')));	

warning('off', 'images:imshow:magnificationMustBeFitForDockedFigure');
dbstop if error
run(fullfile(sun3dsfm_path,'lib/vlfeat/toolbox/vl_setup'));

fprintf('Running SUN3dsfm...\n'); tic;

%% read data
fprintf('Loading data... '); tic;
data = loadConf( base_path, conf_name, frameIDs );
fprintf('Done in %f sec.\n', toc );

if isempty(data)
    fprintf('failed to read file %s\n', conf_name);
    return;
end

%% pairwise transformations
n_matches = length(data.depth) - 1;

fprintf('Computing pairwise transformations...\n');
MatchPairs = cell(1, n_matches);
FrameMatches = zeros( n_matches, 4 );

for frameID = 1:n_matches
  tic;
  next_frameID = min(frameID+1, length(data.depth));
  % try different distances if pair is invalid.
  [MatchPairs{frameID}, output_matches] = align2view(data, frameID, next_frameID, 5.0);
  if(~MatchPairs{frameID}.valid)
    [MatchPairs{frameID}, output_matches] = align2view(data, frameID, next_frameID, 5.5);
    if(~MatchPairs{frameID}.valid)
      [MatchPairs{frameID}, output_matches] = align2view(data, frameID, next_frameID, 6.0);
      if(~MatchPairs{frameID}.valid)
        [MatchPairs{frameID}, output_matches] = align2view(data, frameID, next_frameID, 6.5);
        if(~MatchPairs{frameID}.valid)
          [MatchPairs{frameID}, output_matches] = align2view(data, frameID, next_frameID, 7.0);
          if(~MatchPairs{frameID}.valid)
            [MatchPairs{frameID}, output_matches] = align2view(data, frameID, next_frameID, 7.5);
            if(~MatchPairs{frameID}.valid)
              [MatchPairs{frameID}, output_matches] = align2view(data, frameID, next_frameID, 8.0);
              if(~MatchPairs{frameID}.valid)
                [MatchPairs{frameID}, output_matches] = align2view(data, frameID, next_frameID, 8.5);
              end;
            end;
          end;
        end;
      end;
    end;
  end;
  FrameMatches(frameID, :) = output_matches;
  fprintf('Time to complete frame %d was %f\n', frameID, toc );
  fprintf('----------------------------------------------\n');
end
  

%% naive approach: just put all results together
cameraRtC2W = repmat([eye(3) zeros(3,1)], [1,1,length(data.depth)]);
for frameID = 1:length(data.depth)-1
  cameraRtC2W(:,:,frameID+1) = [cameraRtC2W(:,1:3,frameID) * MatchPairs{frameID}.Rt(:,1:3) cameraRtC2W(:,1:3,frameID) * MatchPairs{frameID}.Rt(:,4) + cameraRtC2W(:,4,frameID)];
end

% save results
matches_folder = fullfile( base_path, 'matches' );
if ~exist(matches_folder, 'dir')
    mkdir( matches_folder );
end

local_matches_name = fullfile( base_path, 'matches/pairwise_matches.txt' )
save_matches( local_matches_name, MatchPairs, data, 1 );
outputConf( fullfile(base_path, strcat(sequence_name, '_pairwise_sift', '.conf')), data, cameraRtC2W );
%outputPly( fullfile(base_path, strcat(sequence_name, '_pairwise_sift.ply')), data, cameraRtC2W );

