function [] = save_matches( file_name, MatchPairs, data, include_all )

  file = fopen( file_name, 'w' );
  n_valid_matches = 0;
  for p = 1:size( MatchPairs, 2 )
    if MatchPairs{ p }.valid || include_all
      n_valid_matches = n_valid_matches+1;
    end
  end
  
  fprintf(file, 'n_matches %d\n', n_valid_matches);
  for p = 1:size( MatchPairs, 2 )
    pair_data = MatchPairs{ p };
    if pair_data.valid || include_all
      matches = pair_data.matches;
      n_matches = size( matches, 2 );
      
      T = [pair_data.Rt; 0 0 0 1];
      C = [ 1 0 0 0 ; 0 -1 0 0 ; 0 0 -1 0 ; 0 0 0 1];
      T = C * T * C;

      Rt = T(1:3, :);

      transformation = reshape( Rt', 1, 12 );

      [~, name_i, ~] = fileparts( data.depth{ pair_data.i } );
      [~, name_j, ~] = fileparts( data.depth{ pair_data.j } );
      scan_name_i = strcat( 'SCAN:', name_i );
      scan_name_j = strcat( 'SCAN:', name_j );
        
      fprintf( file, '%30s %30s %d %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f \n', scan_name_i, scan_name_j, n_matches, transformation );
    end
  end
  fclose( file );
end

