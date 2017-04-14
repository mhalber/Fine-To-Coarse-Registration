// FineToCoarseRegistration - fetregister - io
//
// This file contains functions used for communicating with the storage - reading
// and writing files.

////////////////////////////////////////////////////////////////////////////////
// License:
//
// Copyright (c) 2017 Maciej Halber
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// interface
////////////////////////////////////////////////////////////////////////////////
#ifndef IO_H
#define IO_H

namespace io
{
  i32 read_configuration_file( const char *filename, 
                               SequenceData *seq_data, 
                               Options *opts );

  i32 write_configuration_file( const char *filename, 
                                const SequenceData *seq_data, 
                                const Options *opts );

  i32 write_tum_trajectory_file( const char * filename, 
                                 const SequenceData * seq_data, 
                                 const Options *opts );

  i32 read_sens_file( const char *filename,
                      SequenceData *seq_data, 
                      const Options *opts );

  i32 write_sens_file( const char *filename );

  i32 read_fet_file( const char *filename,
                     SequenceData *seq_data, 
                    Options *opts );

  i32 write_fet_file( const char *filename );

  i32 write_errors( const char * filename,
                     const SequenceData * seq_data );

  i32 write_ply( const char * filename, 
                 const SequenceData * seq_data,
                  const Options * opts,
                  float max_depth = 10.0f );

  i32 write_ply( const char * filename,
                 const FETReconstruction * reconstruction,
                 const Options * opts,
                 float max_depth = 10.0f );
}

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

#ifdef IO_IMPLEMENTATION

// Private Helpers
bool
prv_are_nonnegative( i32 a, i32 b, i32 dummy )
{
  return (a >= 0 && b >= 0);
}

bool
prv_are_divisible( i32 a, i32 b, i32 c)
{
  return ((a % c) == 0 && (b % c) == 0);
}

void
prv_preserve_only_valid_corrs( SequenceData * seq_data, 
                               bool (*is_corr_valid)(i32, i32, i32), 
                               const i32 idx_divisor = 1 )
{
  i32 prev_n_corrs = seq_data->n_gt_corrs;
  i32 new_n_corrs = 0;
  Correspondence * tmp_corrs = (Correspondence*)malloc( prev_n_corrs * 
                                                        sizeof(Correspondence));

  // go through all and only replace ones that are valid, according to validity
  // funct
  for ( i32 corr_idx = 0 ;
        corr_idx < prev_n_corrs;
        corr_idx++ )
  {
    i32 idx_A = seq_data->gt_corrs[ corr_idx ].idx_A;
    i32 idx_B = seq_data->gt_corrs[ corr_idx ].idx_B;

    // check if skip
    if ( !is_corr_valid( idx_A, idx_B, idx_divisor ) )
    {
      continue;
    }

    // copy data
    // if there is discrepancy in terms of shapes in .conf and 
    // .fet we might need to change the indices. Hence 'idx_divisor'
    tmp_corrs[new_n_corrs].name_A  = seq_data->gt_corrs[ corr_idx ].name_A; 
    tmp_corrs[new_n_corrs].name_B  = seq_data->gt_corrs[ corr_idx ].name_B; 
    tmp_corrs[new_n_corrs].idx_A   = seq_data->gt_corrs[ corr_idx ].idx_A / idx_divisor; 
    tmp_corrs[new_n_corrs].idx_B   = seq_data->gt_corrs[ corr_idx ].idx_B / idx_divisor;
    tmp_corrs[new_n_corrs].point_A = seq_data->gt_corrs[ corr_idx ].point_A;
    tmp_corrs[new_n_corrs].point_B = seq_data->gt_corrs[ corr_idx ].point_B;  
    new_n_corrs++;
  }

  if ( new_n_corrs != prev_n_corrs )
  {
    printf("  Warning: Only using subset of gt_correspondences!"
           " ( %d out of %d) \n", new_n_corrs, prev_n_corrs );
  }

  // delete old corrs and copy in from tmp storage
  free( seq_data->gt_corrs );
  seq_data->n_gt_corrs = new_n_corrs;
  i32 copy_size = new_n_corrs * sizeof(Correspondence );
  seq_data->gt_corrs = (Correspondence*)malloc( copy_size );
  memcpy( seq_data->gt_corrs, tmp_corrs, copy_size );

  // free tmp storage
  free(tmp_corrs);

}

////////////////////////////////////////////////////////////////////////////////
// Conf Reading
////////////////////////////////////////////////////////////////////////////////
i32 
parse_intrinsics_cmd( const char *buffer, 
                       char *cmd, 
                       SequenceData *data )
{
  char filename[ MAX_STR_LENGTH ];
  if ( sscanf( buffer, "%s%s", cmd, filename ) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'intrinsics' command format\n" );
    return 0;
  }
  FILE *fp = fopen( filename, "r" );
  if ( !fp )
  {
    fprintf( stderr, "\nUnable to read intrinsics file %s\n", filename );
    return 0;
  }

  bsc::mat3 mat;
  fscanf( fp, "%f %f %f  %f %f %f  %f %f %f",
              &(mat[0][0]), &(mat[1][0]), &(mat[2][0]),
              &(mat[0][1]), &(mat[1][1]), &(mat[2][1]),
              &(mat[0][2]), &(mat[1][2]), &(mat[2][2]) );

  fclose(fp);

  if ( !strcmp( cmd, "color_intrinsics" ) )
  {
    data->color_intrinsics = mat;
    data->color_intrinsics_filename = strdup(filename);
  }
  else if ( !strcmp( cmd, "intrinsics" ) || 
            !strcmp( cmd, "depth_intrinsics" ) ) // NOTE: This should say depth
  {
    data->depth_intrinsics = mat;
    data->depth_intrinsics_filename = strdup(filename);
  }

  return 1;
}

i32 
parse_extrinsics_cmd( const char *buffer, 
                       char *cmd, 
                       SequenceData *data )
{
  char filename[ MAX_STR_LENGTH ];
  if ( sscanf( buffer, "%s%s", cmd, filename ) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'extrinsics' command format\n" );
    return 0;
  }

  FILE *fp = fopen( filename, "r" );
  if ( !fp )
  {
    fprintf( stderr, "\nUnable to read extrinsics file %s\n", filename );
    return 0;
  }

  bsc::mat4 mat;
  fscanf( fp, "%f %f %f %f  %f %f %f %f  %f %f %f %f  %f %f %f %f",
              &(mat[0][0]), &(mat[1][0]), &(mat[2][0]), &(mat[3][0]),
              &(mat[0][1]), &(mat[1][1]), &(mat[2][1]), &(mat[3][1]),
              &(mat[0][2]), &(mat[1][2]), &(mat[2][2]), &(mat[3][2]),
              &(mat[0][3]), &(mat[1][3]), &(mat[2][3]), &(mat[3][3]) );

  fclose(fp);

  if ( !strcmp( cmd, "color_extrinsics" ) )
  {
    data->color_extrinsics = mat;
    data->color_extrinsics_filename = strdup(filename);
  }
  else if ( !strcmp( cmd, "depth_extrinsics" ) )
  {
    data->depth_extrinsics = mat;
    data->depth_extrinsics_filename = strdup(filename);
  }
  
  return 1;
}

i32 
parse_extrinsics_rgbdsfm_cmd( const char *buffer, 
                           char *cmd, 
                           SequenceData *data )
{
  char filename[ MAX_STR_LENGTH ];
  if ( sscanf( buffer, "%s%s", cmd, filename ) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'extrinsics_rgbdsfm' command format\n" );
    return 0;
  }

  FILE *fp = fopen( filename, "r" );
  if ( !fp )
  {
    fprintf( stderr, "\nUnable to read rgbdsfm extrinsics file %s\n", 
             filename );
    return 0;
  }
  else {
    bsc::mat4d pre_align = bsc::scale(bsc::mat4d(), bsc::vec3d( -1.0, 1.0, -1.0 ) );
    data->rgbdsfm_xforms = (bsc::mat4d*)
                                   malloc( sizeof(data->rgbdsfm_xforms[0]) * data->n_frames );
    for ( i32 idx = 0 ; idx < data->n_frames ; ++idx )
    {
      bsc::mat4d mat;
      fscanf( fp, "%lf %lf %lf %lf  %lf %lf %lf %lf  %lf %lf %lf %lf",
                &(mat[0][0]), &(mat[1][0]), &(mat[2][0]), &(mat[3][0]),
                &(mat[0][1]), &(mat[1][1]), &(mat[2][1]), &(mat[3][1]),
                &(mat[0][2]), &(mat[1][2]), &(mat[2][2]), &(mat[3][2]) );
      mat = bsc::scale( mat, bsc::vec3d( 1.0, -1.0, -1.0 ) );
      data->rgbdsfm_xforms[idx] = pre_align * mat;
    }
    fclose(fp);
  }
  return 1;
}

i32 
parse_extrinsics_robust_reconstruction_cmd( const char *buffer, 
                           char *cmd, 
                           SequenceData *data,
                           Options *opts )
{
  char filename[ MAX_STR_LENGTH ];
  if ( sscanf( buffer, "%s%s", cmd, filename ) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'extrinsics_robust_reconstruction' "
                     "command format\n" );
    return 0;
  }

  FILE *fp = fopen( filename, "r" );
  if ( !fp )
  {
    fprintf( stderr, "\nUnable to read robust reconstruction "
                     "extrinsics file %s\n", filename );
    return 0;
  }
  else {
    bsc::mat4d pre_align = bsc::scale(bsc::mat4d(), bsc::vec3d( -1.0, 1.0, -1.0 ) );
    data->robust_reconstruction_xforms = (bsc::mat4d*) 
                                   malloc( sizeof(data->robust_reconstruction_xforms[0]) * data->n_frames );
    // read in existing frames
    i32 last_frame = -1;
    for ( i32 idx = 0 ; idx < data->n_frames ; ++idx )
    {
      i32 id1, id2, frame;
      bsc::mat4d mat;
      fscanf(fp, "%d %d %d", &id1, &id2, &frame);
      // robust reconstruction throws some frames out
      if ( idx > frame )
      {
        break;
      }
      fscanf(fp, "%lf %lf %lf %lf", &mat[0][0], &mat[1][0], &mat[2][0], &mat[3][0]);
      fscanf(fp, "%lf %lf %lf %lf", &mat[0][1], &mat[1][1], &mat[2][1], &mat[3][1]);
      fscanf(fp, "%lf %lf %lf %lf", &mat[0][2], &mat[1][2], &mat[2][2], &mat[3][2]);
      fscanf(fp, "%lf %lf %lf %lf", &mat[0][3], &mat[1][3], &mat[2][3], &mat[3][3]);
      mat = bsc::scale( mat, bsc::vec3d( 1.0, -1.0, -1.0 ) );

      data->robust_reconstruction_xforms[idx] = mat;
      last_frame = frame;
    }
    // make the rest zero
    for ( i32 idx = last_frame ; idx < data->n_frames ; ++idx )
    {
      data->robust_reconstruction_xforms[idx] = pre_align * bsc::mat4d(0.0f);
    }
    fclose(fp);
  }
  return 1;
}

i32 
parse_extrinsics_tum_format_cmd( const char *buffer, 
                                  char *cmd, 
                                  SequenceData *data,
                                  bsc::mat4d ** xforms_ptr )
{
  char filename[ MAX_STR_LENGTH ];
  if ( sscanf( buffer, "%s%s", cmd, filename ) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'extrinsics elastic_fusion/kintinuous'"
                     " command format\n" );
    return 0;
  }

  FILE *fp = fopen( filename, "r" );
  if ( !fp )
  {
    fprintf( stderr, "\nUnable to read elastic fusion/kintinuous "
                     "extrinsics file %s\n", filename );
    return 0;
  }
  else {
    bsc::mat4d pre_align = bsc::scale(bsc::mat4d(), bsc::vec3d( -1.0, 1.0, -1.0 ) );
    bsc::vec3d flip = bsc::vec3d(1.0f, 1.0f, 1.0f);
    // Used for kintinuous and elastic fusion for sun3d scens
    if ( strcmp(data->dataset, "sun3d") == 0 )
    {
      flip = bsc::vec3d( 1.0f, -1.0f, -1.0f );
    }
    else if( strcmp( data->dataset, "icl") == 0 )
    {
      flip = bsc::vec3d( 1.0f, 1.0f, -1.0f );
    }
        else if( strcmp( data->dataset, "tum") == 0 )
    {
      flip = bsc::vec3d( 1.0f, 1.0f, -1.0f );
    }
    
    (*xforms_ptr) = (bsc::mat4d*) malloc( sizeof(bsc::mat4d) * data->n_frames );
    for ( i32 idx = 0 ; idx < data->n_frames ; ++idx )
    {
      bsc::quatd orientation;
      bsc::vec3d position;
      r32 timestamp;

      if ( fscanf( fp, "%f  %lf %lf %lf  %lf %lf %lf %lf\n",
                       &timestamp,
                       &position.x,
                       &position.y,
                       &position.z,
                       &orientation.x,
                       &orientation.y,
                       &orientation.z,
                       &orientation.w ) != 8 ) 
      {
         // not sure what to do...
      }
      bsc::mat4d mat = bsc::to_mat4( orientation );
      mat[3] = bsc::vec4d( position, 1 );
      mat = bsc::scale( mat, flip );
      (*xforms_ptr)[idx] = pre_align * mat;
    }
    fclose(fp);
  }
  return 1;
}


i32
parse_pairwise_matches_cmd( const char *buffer, 
                            char *cmd, 
                            SequenceData *data )
{

  char filename[MAX_STR_LENGTH];
  if( sscanf( buffer, "%s%s", cmd, filename) != (unsigned int) 2 )
  {
    fprintf( stderr, "Invalid 'pairwise_matches' command format\n" );
    return 0;
  }
  data->pairwise_matches_filename = strdup( filename );
  // int sift_corrs_threshold = 10;

  FILE * file = fopen( filename, "r" );
  if( file )
  {
    size_t n_matches = 0;
    fscanf( file, "n_matches %lu\n", &n_matches );
    if( n_matches + 1 != (size_t)data->n_frames ) 
    {
      fprintf(stderr, "Number of pairwise matches is incorrect (%lu vs. %d)!\n", n_matches + 1, data->n_frames );
    }
    data->pairwise_xforms = (bsc::mat4d*)malloc( data->n_frames * 
                                              sizeof(data->pairwise_xforms[0]));
    data->pairwise_xforms_inliers = (int *)malloc( data->n_frames * 
                                                   sizeof(int));
    // set all as identity initially
    for ( i32 i = 0; i < data->n_frames ; i++ )
    {
      data->pairwise_xforms[i] = bsc::mat4d();
      data->pairwise_xforms_inliers[i] = 0;
    }

    int min_corrs = 1e9;
    int max_corrs = -1e9;
    float sum_of_corrs = 0;
    for ( i32 i = 1 ; 
          i < data->n_frames ;
          ++i )
    {
      i32 n_sift_corrs;
      char name_A[ MAX_STR_LENGTH ];
      char name_B[ MAX_STR_LENGTH ];
      bsc::mat4d m;
      fscanf( file, "%s %s %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &name_A[0], &name_B[0], &n_sift_corrs,
                    &m[0][0], &m[1][0], &m[2][0], &m[3][0],
                    &m[0][1], &m[1][1], &m[2][1], &m[3][1],
                    &m[0][2], &m[1][2], &m[2][2], &m[3][2] );
      
      // if( bsc::norm(bsc::vec3d(m[3])) > 0.05 ) { m = bsc::mat4d(); n_sift_corrs = 10.0; } /* */

      data->pairwise_xforms[i] = data->pairwise_xforms[i-1] * m;
      data->pairwise_xforms_inliers[i-1] = n_sift_corrs;

      if ( n_sift_corrs > max_corrs ) max_corrs = n_sift_corrs;
      if ( n_sift_corrs < min_corrs ) min_corrs = n_sift_corrs;
      sum_of_corrs += (float) n_sift_corrs;
    }
    data->pairwise_xforms_inliers[data->n_frames-1] = data->pairwise_xforms_inliers[data->n_frames-2];
     
    fclose( file );
  }
    else
  {
    printf("Could not open pairwise matches file %s!\n", filename );
    return 0;
  }

  return 1;
}



i32 
parse_n_correspondences_cmd( const char *buffer, char *cmd, SequenceData *data )
{
  i32 n_gt_corrs = 0;
  if ( sscanf( buffer, "%s%d", cmd, &n_gt_corrs) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'n_correspondences' command format\n" );
    return 0;
  }

  data->n_gt_corrs = n_gt_corrs;

  if ( data->gt_corrs== NULL )
  {
    data->gt_corrs= (Correspondence*)malloc( n_gt_corrs * sizeof(Correspondence) );
  }
  return 1;
}

i32 
parse_point3d_cmd( const char *buffer, 
                   char *cmd, 
                   SequenceData *data,
                   i32 idx )
{
  char name_A[MAX_STR_LENGTH], name_B[MAX_STR_LENGTH];
  bsc::vec3d p_A, p_B;
  if ( sscanf( buffer, "%s %s %s %lf %lf %lf %lf %lf %lf", cmd, 
                       name_A, name_B, 
                       &(p_A.x), &(p_A.y), &(p_A.z),
                       &(p_B.x), &(p_B.y), &(p_B.z) ) != (u32) 9 )
  {
    fprintf( stderr, "Invalid 'point3d' command format\n" );
    return 0;
  }

  Correspondence * c = &(data->gt_corrs[ idx ]);
  c->name_A = strdup( name_A );
  c->name_B = strdup( name_B );
  c->point_A = p_A;
  c->point_B = p_B; 

  return 1;
}

i32
parse_correspondences_cmd( const char *buffer, 
                           char *cmd, 
                           SequenceData *seq_data )
{

  char filename[MAX_STR_LENGTH];
  if ( sscanf( buffer, "%s%s", cmd, filename) != (unsigned int) 2 )
  {
    fprintf( stderr, "Invalid 'correspondences' command format\n" );
    return 0;
  }
  seq_data->correspondences_filename = strdup( filename );

  FILE * file = fopen( filename, "r" );
  if ( file )
  {
    char buffer[MAX_STR_LENGTH];
    bool success = 1;
    i32 n_read_corrs = 0;
    while ( fgets( buffer, MAX_STR_LENGTH, file ) )
    {
      char cmd[MAX_STR_LENGTH];
      if ( sscanf(buffer, "%s", cmd) != ( u32) 1)
        continue;

      if (cmd[0] == '#')
        continue;

      if ( !strcmp(cmd, "n_correspondences" ) )
      {
        success &= parse_n_correspondences_cmd( buffer, cmd, seq_data );
      }
      else if ( !strcmp(cmd, "point3d" ) )
      {
        success &= parse_point3d_cmd( buffer, cmd, seq_data, n_read_corrs++ );
      }
      else
      {
        continue;
      }
    }

    if ( !success )
    {
      return 0;
    }
  }
  else
  {
    printf("Could not open correspondence file %s!\n", filename );
    return 0;
  }

  return 1;
}

i32 
parse_dataset_cmd( const char *buffer, 
                    char *cmd, 
                    SequenceData *data )
{
  char dataset[256];
  if ( sscanf( buffer, "%s%s", cmd, dataset ) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'dataset' command format\n" );
    return 0;
  }
  data->dataset = strdup( dataset );
  return 1;
}


i32 
parse_resolution_cmd( const char *buffer, 
                             char *cmd, 
                             SequenceData *data )
{
  i32 x, y;
  if ( sscanf( buffer, "%s%d%d", cmd, &x, &y ) != (u32) 3 )
  {
    fprintf( stderr, "Invalid 'resolution' command format\n" );
    return 0;
  }

  if ( !strcmp( cmd, "color_resolution" ) )
  {
    data->color_res = bsc::vec2i( x, y );
  }
  else if ( !strcmp( cmd, "depth_resolution" ) )
  {
    data->depth_res = bsc::vec2i( x, y );
  }
  return 1;
}


i32
parse_dir_cmd ( const char * buffer, char * cmd, SequenceData * data )
{
  char dirname[MAX_STR_LENGTH];

  if ( sscanf( buffer, "%s%s", cmd, dirname) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'directory' command format\n" );
    return 0;
  } 
  if ( strcmp( cmd, "color_directory" ) == 0 ||
       strcmp( cmd, "image_directory") == 0 )
  {
    data->color_dir = strdup( dirname );
  }
  else if ( !strcmp( cmd, "depth_directory" ) )
  {
    data->depth_dir = strdup( dirname );
  }
  return 1;
}

i32 
parse_n_frames_cmd( const char *buffer, char *cmd, SequenceData *data )
{
  i32 n_images = 0;
  if ( sscanf( buffer, "%s%d", cmd, &n_images) != (u32) 2 )
  {
    fprintf( stderr, "Invalid 'n_frames' command format\n" );
    return 0;
  }

  data->n_frames = n_images;

  return 1;
}

i32 
parse_scan_cmd( const char *buffer, char *cmd, SequenceData *data, i32 & idx )
{
  if ( idx >= data->n_frames ) return 1;
  char depth_name[256], color_name[256];
  bsc::mat4d mat;
  if (sscanf( buffer, "%s %s %s"
                      "%lf %lf %lf %lf"
                      "%lf %lf %lf %lf"
                      "%lf %lf %lf %lf"
                      "%lf %lf %lf %lf", 
          cmd, depth_name, color_name,
          &(mat[0][0]), &(mat[1][0]), &(mat[2][0]), &(mat[3][0]),
          &(mat[0][1]), &(mat[1][1]), &(mat[2][1]), &(mat[3][1]),
          &(mat[0][2]), &(mat[1][2]), &(mat[2][2]), &(mat[3][2]),
          &(mat[0][3]), &(mat[1][3]), &(mat[2][3]), &(mat[3][3]) ) != (u32) 19 )
  {
    fprintf( stderr, "Invalid %d-th 'scan' command format\n", idx );
    return 0;
  }

  if ( data->depth_names == NULL )
  {
    data->depth_names = (char(*)[MAX_FRAME_NAME_LENGTH])malloc( data->n_frames * MAX_FRAME_NAME_LENGTH );
  }
  if ( data->color_names == NULL )
  {
    data->color_names = (char(*)[MAX_FRAME_NAME_LENGTH])malloc( data->n_frames * MAX_FRAME_NAME_LENGTH );
  }
  if ( data->xforms_history == NULL )
  {
    data->xforms_history = (bsc::mat4d*)malloc(data->n_frames * MAX_ITER * sizeof( data->xforms_history[0] ) );
    data->n_xforms = 1;
  }

  strcpy( data->color_names[idx], color_name );
  strcpy( data->depth_names[idx], depth_name );
  data->xforms_history[ idx ] = mat;
  idx++;

  return 1;
}


i32 io::
read_configuration_file( const char * filename, 
                         SequenceData * seq_data, 
                         Options *opts )
{
  bsc::stop_watch timer;
  timer.read();

  const char *conf_filename = filename;

  FILE *configuration_fp = fopen( conf_filename, "r" );
  if ( !configuration_fp )
  {
    fprintf(stderr, "Unable to open configuration file %s\n",
                     conf_filename);
    exit(-1);
  }

  if ( opts->print_verbose )
  {
    printf("Reading configuration file...\n");
  }

   // Parse file
  char buffer[MAX_STR_LENGTH];
  i32 line_number = 0;

  bool success = 1;
  i32 n_read_frames = 0;
  while ( fgets( buffer, MAX_STR_LENGTH, configuration_fp ) )
  {
    char cmd[MAX_STR_LENGTH];
    line_number++;
    if ( sscanf(buffer, "%s", cmd) != ( u32) 1)
      continue;

    if (cmd[0] == '#')
      continue;

    if ( !strcmp(cmd, "dataset" ) )
    {
      success &= parse_dataset_cmd( buffer, cmd, seq_data );
    }
    else if ( !strcmp(cmd, "intrinsics") )
    {
      success &= parse_intrinsics_cmd( buffer, cmd, seq_data );
    }
    else if ( !strcmp(cmd, "depth_intrinsics") ) 
    {
      success &= parse_intrinsics_cmd( buffer, cmd, seq_data );
    }
    else if ( !strcmp(cmd, "color_intrinsics" ) )
    {
      success &= parse_intrinsics_cmd( buffer, cmd, seq_data );
    }
    else if ( !strcmp(cmd, "depth_extrinsics" ) )
    {
      success &= parse_extrinsics_cmd( buffer, cmd, seq_data );
    }
    else if ( !strcmp(cmd, "color_extrinsics" ) )
    {
      success &= parse_extrinsics_cmd(buffer, cmd, seq_data );
    }
    else if ( !strcmp(cmd, "extrinsics_rgbdsfm" ) )
    {
      // OPTIONAL
      parse_extrinsics_rgbdsfm_cmd(buffer, cmd, seq_data );
    }
    else if ( !strcmp(cmd, "extrinsics_robust_reconstruction" ) )
    {
      // OPTIONAL 
      parse_extrinsics_robust_reconstruction_cmd(buffer, cmd, seq_data, opts );
    }
    else if ( !strcmp(cmd, "extrinsics_elastic_fusion" ) )
    {
      //OPTIONAL
      parse_extrinsics_tum_format_cmd(buffer, cmd, seq_data, 
                                           &(seq_data->elastic_fusion_xforms) );
    }
    else if ( !strcmp(cmd, "extrinsics_kintinuous" ) )
    {
      //OPTIONAL
      parse_extrinsics_tum_format_cmd(buffer, cmd, seq_data,
                                               &(seq_data->kintinuous_xforms) );
    }
    else if ( !strcmp(cmd, "extrinsics_gt" ) )
    {
      // TODO(maciej): Add switch based on the dataset
      parse_extrinsics_tum_format_cmd(buffer, cmd, seq_data,
                                               &(seq_data->gt_xforms) );
    }
    else if ( !strcmp(cmd, "depth_resolution" ) )
    {
      success &= parse_resolution_cmd( buffer, cmd, seq_data );
    }
    else if ( !strcmp(cmd, "color_resolution" ) )
    {
      success &= parse_resolution_cmd(buffer, cmd, seq_data );
    }
    else if ( !strcmp( cmd, "depth_directory" ) ) 
    {
      success &= parse_dir_cmd( buffer, cmd, seq_data );
    }
    else if ( strcmp( cmd, "image_directory" ) == 0 ) 
    {
      success &= parse_dir_cmd( buffer, cmd, seq_data );
    }
    else if ( strcmp( cmd, "color_directory" ) == 0 ) 
    {
      success &= parse_dir_cmd( buffer, cmd, seq_data );
    }
    else if ( strcmp( cmd, "n_images" ) == 0 || strcmp( cmd, "n_frames" ) == 0   )
    {
      success &= parse_n_frames_cmd(buffer, cmd, seq_data );
    }
    else if (!strcmp( cmd, "pairwise_matches" ) )
    {
      success &= parse_pairwise_matches_cmd(buffer, cmd, seq_data );
    }
    else if (!strcmp( cmd, "correspondences" ) )
    {
      success &= parse_correspondences_cmd(buffer, cmd, seq_data );
    }
    else if (!strcmp(cmd, "scan"))
    {
      success &= parse_scan_cmd( buffer, cmd, seq_data, n_read_frames );
    }
  }

  if ( !success )
  {
    return 0;
  }

  // set the transformation pointers in place
  seq_data->current_xforms = seq_data->xforms_history;
  seq_data->initial_xforms = seq_data->xforms_history;
  seq_data->format = "conf";
  if ( !seq_data->dataset ) seq_data->dataset = (char*)"Unknown";

  // override missing / incorrect options
  if ( opts->frame_end == 0 ||
       opts->frame_end > seq_data->n_frames  ) 
  { 
    opts->frame_end = seq_data->n_frames; 
  }

  if ( seq_data->n_frames > n_read_frames )
  {
    opts->frame_end = n_read_frames;
  }

  if ( seq_data->color_res == bsc::vec2i(0, 0) )
  {
    seq_data->color_res = bsc::vec2i( 640, 480 );
  }

  if ( seq_data->depth_res == bsc::vec2i(0, 0) )
  {
    seq_data->depth_res = bsc::vec2i( 640, 480 );
  }

  // find frame indices for correspondences
  // I dislike this greatly
  // Some hash map should replace this
  int n_valid = 0;
  for ( int corr_idx = 0; corr_idx < seq_data->n_gt_corrs ; corr_idx++ )
  {
    char * name_A = seq_data->gt_corrs[ corr_idx ].name_A;
    char * name_B = seq_data->gt_corrs[ corr_idx ].name_B;
    seq_data->gt_corrs[ corr_idx ].idx_A = -1;
    seq_data->gt_corrs[ corr_idx ].idx_B = -1;
    for ( i32 frame_idx = 0;
          frame_idx < seq_data->n_frames;
          frame_idx++ )
    {
      char scan_name[MAX_STR_LENGTH];
      sprintf( scan_name, "SCAN:%s", seq_data->depth_names[ frame_idx ] );
      i32 len = strlen( scan_name );
      scan_name[len-4] = '\0';
      if ( !strcmp( scan_name, name_A ) )
      {
        seq_data->gt_corrs[ corr_idx ].idx_A = frame_idx;
      }
      if ( !strcmp( scan_name, name_B ) )
      {
        seq_data->gt_corrs[ corr_idx ].idx_B = frame_idx;
      }
    }

    if ( seq_data->gt_corrs[ corr_idx ].idx_A != -1 &&
         seq_data->gt_corrs[ corr_idx ].idx_B != -1 )
    {
      n_valid++;
    }

  }
  // if there are some invalid corrs, lets remove them
  if ( n_valid != seq_data->n_gt_corrs )
  {
    prv_preserve_only_valid_corrs( seq_data, 
                                   prv_are_nonnegative );
  }

  // Close configuration file
  fclose(configuration_fp);

  // Print statistics
  if ( opts->print_verbose )
  {
    printf("Done in %5.3f sec.\n", timer.elapsed() );
    printf("  # Frames = %d(%d)\n", seq_data->n_frames, n_read_frames );
    printf("  # Correspondences = %d\n", seq_data->n_gt_corrs );
    printf("  # Color Directory = %s\n", seq_data->color_dir );
    printf("  # Depth Directory = %s\n", seq_data->depth_dir );
    printf("\n");
  }

  // Return success
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
// Conf Writing
////////////////////////////////////////////////////////////////////////////////

i32 io::
write_configuration_file( const char * filename, 
                          const SequenceData * seq_data, 
                          const Options *opts )
{
  FILE * fp = fopen( filename, "w" );
  if ( fp )
  {
    fprintf( fp, "dataset %s\n", seq_data->dataset );
    fprintf( fp, "intrinsics %s\n", seq_data->depth_intrinsics_filename );
    if ( seq_data->color_intrinsics_filename )
    {
      fprintf( fp, "color_intrinsics %s\n", seq_data->color_intrinsics_filename );
    }
    fprintf( fp, "n_images %d\n", seq_data->n_frames );
    fprintf( fp, "depth_directory %s\n", seq_data->depth_dir );
    fprintf( fp, "image_directory %s\n", seq_data->color_dir );
    fprintf( fp, "color_resolution %d %d\n", seq_data->color_res.x, seq_data->color_res.y );
    fprintf( fp, "depth_resolution %d %d\n", seq_data->depth_res.x, seq_data->depth_res.y );
    if ( seq_data->correspondences_filename )
    {
      fprintf( fp, "correspondences %s\n", seq_data->correspondences_filename );
    }
    if ( seq_data->pairwise_matches_filename )
    {
      fprintf( fp, "pairwise_matches %s\n", seq_data->pairwise_matches_filename );
    }

    fprintf( fp, "extrinsics_rgbdsfm extrinsics/rgbdsfm.txt\n" );
    fprintf( fp, "extrinsics_robust_reconstruction extrinsics/robust_reconstruction.txt\n" );
    fprintf( fp, "extrinsics_elastic_fusion extrinsics/elastic_fusion_poses_improved.txt\n" );
    fprintf( fp, "extrinsics_kintinuous extrinsics/kintinuous_poses.txt\n" );
    
    fprintf( fp, "\n" );

    for ( i32 i = 0 ; i < opts->frame_end; i++ )
    {
      bsc::mat4d mat = seq_data->current_xforms[ i ];
      fprintf( fp, "scan %30s %30s  %lf %lf %lf %lf  %lf %lf %lf %lf  %lf %lf %lf %lf  %lf %lf %lf %lf\n",
                    seq_data->depth_names[i], seq_data->color_names[i],
                    mat[0][0], mat[1][0], mat[2][0], mat[3][0],
                    mat[0][1], mat[1][1], mat[2][1], mat[3][1],
                    mat[0][2], mat[1][2], mat[2][2], mat[3][2],
                    mat[0][3], mat[1][3], mat[2][3], mat[3][3] ); 
    }
    fclose( fp );
  }
  else
  {
    printf( "Could not open file %s\n", filename );
    return 0;
  }
  return 1;
}


i32 io::
write_tum_trajectory_file( const char * filename, 
                           const SequenceData * seq_data, 
                           const Options *opts )
{
  FILE * fp = fopen( filename, "w" );
  if ( fp )
  {
    for ( i32 i = 0 ; i < opts->frame_end; i++ )
    {
      std::string name ( seq_data->depth_names[i] );
      size_t last_index = name.find_last_of("."); 
      std::string raw_name = name.substr(0, last_index); 
    r64 depth_timestamp = atof( raw_name.c_str() );

      bsc::mat4d mat = seq_data->current_xforms[ i ];
      bsc::vec3d pos = bsc::vec3d( mat[3] );
      bsc::quatd q = bsc::to_quat( mat );

      fprintf( fp, "%f  %lf %lf %lf  %lf %lf %lf %lf\n",
                    depth_timestamp,
                    pos.x, pos.y, pos.z,
                    q.x, q.y, q.z, q.w ); 
    }
    fclose( fp );
  }
  else
  {
    printf( "Could not open file %s\n", filename );
    return 0;
  }
  return 1;
}

////////////////////////////////////////////////////////////////////////////////
// FET reading
////////////////////////////////////////////////////////////////////////////////

i32 io::
read_fet_file( const char * filename,
             SequenceData * seq_data, 
             Options *opts )
{
  bsc::stop_watch timer;
  timer.read();
  if ( opts->print_verbose )
  {
    printf("Reading reconstruction file...\n");
  }

  // The whole FET format is a bit meh.
  FETReconstruction * tmp_reconstruction       = new FETReconstruction();
  seq_data->reconstruction = new FETReconstruction();
  if ( !tmp_reconstruction->ReadFile( filename ) )
  {
    printf( "Could not open %s file!\n", filename );
    return 0;
  }

  // Set solver
  if (opts->solver >= 0) seq_data->reconstruction->solver = opts->solver;

  // Let's store some stuff
  seq_data->format   = "FET";
  if ( seq_data->dataset == NULL )
  {
    seq_data->dataset = (char*)"Not Stored!";
  }

  // Move shapes from the tmp_reconstruction to actual one
  int i = 0, j = 0;
  char cur_name[512];
  for ( ; i < opts->frame_end ; i++ )
  {
    char * depth_name = seq_data->depth_names[ i ];
    sprintf( cur_name, "SCAN:%s", depth_name ); // add prefix
    char * pch = strchr( cur_name, '.' );       // find '.'
    cur_name[ pch-(&cur_name[0]) ] = '\0';      // terminate at '.' ( remove extension )
    
    FETShape * shape = tmp_reconstruction->Shape(j);
    const char * shape_name = shape ? shape->Name() : NULL;

    if ( shape_name && !strcmp( shape_name, cur_name ) )
    {
      tmp_reconstruction->RemoveShape( shape );
      seq_data->reconstruction->InsertShape( shape );
      j += 1;
    }
    else
    {
      FETShape * shape = new FETShape( NULL );
      shape->SetName( cur_name );
      shape->SetOrigin( R3zero_point );
      seq_data->reconstruction->InsertShape( shape );
    }
  }

  delete tmp_reconstruction;

/*
  // Fill in matches
  for ( int i = 0 ; i < seq_data->n_loop_closure_matches ; ++i )
  {
    Match * m = &seq_data->loop_closure_matches[i];
    FETShape * shape_A = seq_data->reconstruction->Shape(m->name_A);
    FETShape * shape_B = seq_data->reconstruction->Shape(m->name_B);
    if ( shape_A && shape_B )
    {
      m->idx_A = shape_A->reconstruction_index;
      m->idx_B = shape_B->reconstruction_index;
    }
    else
    {
      m->idx_A = -1;
      m->idx_B = -1;
    }
  }
*/

  // report status
  if ( opts->print_verbose )
  {
    printf("Done in %5.3f sec.\n", timer.elapsed() );
    printf("  # Shapes = %d\n", seq_data->reconstruction->NShapes() );
    printf("  # Features = %d\n", seq_data->reconstruction->NFeatures() );
    // printf("  Pose Swap: %d\n", swap_poses );
    printf("\n");
  }

  return 1;
}

i32 io::
read_sens_file( const char *filename,
               SequenceData *seq_data, 
                const Options *opts )
{
  return 0;
}

i32 io::
write_sens_file( const char * filename )
{
  return 0;
}

i32 io::
write_errors( const char *filename,
             const SequenceData * seq_data )
{
  FILE *error_f = fopen( filename, "w" );
  if ( !error_f )
  {
    bsc::error( "IO Error",__LINE__, "Could not open error file!" );
    return 0;
  }
  else
  {
    const char* method_names[5] = { "Ours", "Sun3D", "Robust_Reconstruction", 
                                    "Elastic_Fusion", "Kintinuous" };
    r32 rmse[5], mean[5], sdev[5];
    eval::calculate_error_stats( seq_data, &(rmse[0]), &(mean[0]), &(sdev[0]) );

    for ( int i = 0 ; i < 5 ; ++i )
    {
      fprintf( error_f, "%s %f %f %f\n", method_names[i], rmse[i], mean[i], sdev[i] );
    }
  }

  for ( int i = 0 ; i < seq_data->n_xforms ; i++ )
  {
    fprintf(error_f, "%f, ", seq_data->errors[i] );
  }
  fprintf(error_f, "\n" );
  fclose( error_f );
  return 1;
}

#ifdef FETREGISTER_USE_WINDOW
void 
prv_init_pointcloud( bsc::geometry_data * pointcloud,
                     const bsc::img *color_img,
                     const bsc::img_u16 *depth_img,
                     const bsc::mat4d T,
                     const bsc::mat3 K,
                     const i32 downsample_factor,
                     const bool bitshift,
                     const r32 near_dist,
                     const r32 far_dist )
{
  i32 valid_i = 0;
  i32 w = depth_img->width;
  i32 h = depth_img->height;
  r32 ratio_x = (r32)w / color_img->width; 
  r32 ratio_y = (r32)h / color_img->height; 

  // Do bitshift if needed
  bsc::img_u16 depth_cpy( *depth_img );
  if (bitshift)
  {
    for ( int i = 0 ; i < w * h ; ++ i )
    {
      u16 d = *(depth_img->at(i));
      (*depth_cpy.at(i)) =  ((d >> 3) & 0x1FFF) | ((d & 0x7) << 13);
    }
  }

  // create r2 grid ( TEMP! )
  R2Grid grd( w, h );
  for ( int i = 0 ; i < w * h ; ++ i )
  {
    grd.SetGridValue( i, *(depth_cpy.at(i)) * 0.001f );
  }

  // find pixels at depth discontinuities
  for (i32 y = 0; y < h; y++ )
  {
    for (i32 x = 0; x < w ; x++ )
    {
      if ( x == 0 || x == w-1 ||
           y == 0 || y == h - 1 )
      {
        grd.SetGridValue(x,y,R2_GRID_UNKNOWN_VALUE);
        continue;
      }
      RNScalar d0 = grd.GridValue(x,y);
      RNScalar d1 = grd.GridValue(x+1,y);
      RNScalar d2 = grd.GridValue(x,y+1);
      if ( d0 < 0.0001 || 
           fabs(d1-d0) > 0.1 || 
           fabs(d2-d0) > 0.1 )
      {
        grd.SetGridValue(x,y,R2_GRID_UNKNOWN_VALUE);
      }
    }
  }

  // Erode
  R2Grid mask( grd );
  mask.Substitute(R2_GRID_UNKNOWN_VALUE, 0.0);
  mask.Erode(3.0);

  // Apply errosion mask
  grd.Multiply( mask );
  grd.Substitute(0.0, R2_GRID_UNKNOWN_VALUE);
  
  // Downsample
  grd.Resample( w / downsample_factor, h / downsample_factor );
  grd.Substitute( R2_GRID_UNKNOWN_VALUE, 0.0 );

  w = grd.XResolution();
  h = grd.YResolution();
  ratio_x = color_img->width / (r32)w; 
  ratio_y = color_img->height / (r32)h; 

  bsc::mat3 L( K );
  L[2][0] = (1.0f / ratio_x) * K[2][0];
  L[2][1] = (1.0f / ratio_y) * K[2][1];
  L[0][0] = (1.0f / ratio_x) * K[0][0];
  L[1][1] = (1.0f / ratio_y) * K[1][1];

  bsc::img_u16 depth_final( w, h, 1 );
  for ( int i = 0 ; i < w * h ; ++i )
  {
    (*depth_final.at(i)) = grd.GridValue(i) * 1000.0f;
  }

  for (i32 y = 0; y < h-1; y++)
  {
    for (i32 x = 0; x < w-1; x++)
    {
      const u8 *color = color_img->at( ratio_x * x, ratio_y * y );
      u16 d0 = *(depth_final.at(x, y));
      u16 d1 = *(depth_final.at(x+1, y));
      u16 d2 = *(depth_final.at(x, y+1));
      
      r32 depth0 = d0 * 0.001f;
      r32 depth1 = d1 * 0.001f;
      r32 depth2 = d2 * 0.001f;

      if (depth0 < near_dist || depth0 > far_dist ||
          depth1 < near_dist || depth1 > far_dist ||
          depth2 < near_dist || depth2 > far_dist)
      {
        continue;
      }

      // avoid discontinuities
      if ( fabs( depth0 - depth1 ) > 0.1f ) continue;
      if ( fabs( depth0 - depth2 ) > 0.1f ) continue;

      bsc::vec3d p0( ((x) - L[2][0]) * depth0 / L[0][0],
                    ((h - y) - L[2][1]) * depth0 / L[1][1],
                    -depth0 );
      bsc::vec3d p1( ((x + 1) - L[2][0]) * depth1 / L[0][0],
                    ((h - y) - L[2][1]) * depth1 / L[1][1],
                     -depth1 );
      bsc::vec3d p2( ((x) - L[2][0]) * depth2 / L[0][0],
                    ((h - y - 1) - L[2][1]) * depth2 / L[1][1],
                    -depth2 );

      bsc::vec3d v0 = bsc::normalize( p1 - p0 );
      bsc::vec3d v1 = bsc::normalize( p2 - p0 );

      bsc::vec3d n = bsc::normalize( bsc::cross(v0, v1) );
      bsc::vec3d ray = bsc::normalize( bsc::vec3d( 0.0f, 0.0f, 0.0f ) - p0 );
      if ( bsc::dot(n, ray) < 0 ) n = -n; 

      p0 = bsc::vec3d( T * bsc::vec4d(p0, 1.0) );
      n  = bsc::vec3d( T * bsc::vec4d(n, 0.0) );

      pointcloud->positions[3 * valid_i + 0] = p0.x;
      pointcloud->positions[3 * valid_i + 1] = p0.y;
      pointcloud->positions[3 * valid_i + 2] = p0.z;

      pointcloud->normals[3 * valid_i + 0] = n.x;
      pointcloud->normals[3 * valid_i + 1] = n.y;
      pointcloud->normals[3 * valid_i + 2] = n.z;

      pointcloud->colors_a[3 * valid_i + 0] = color[0];
      pointcloud->colors_a[3 * valid_i + 1] = color[1];
      pointcloud->colors_a[3 * valid_i + 2] = color[2];

      valid_i++;
    }
  }

  pointcloud->n_vertices = valid_i;
}
#endif


i32 io::
write_ply( const char * filename, 
          const FETReconstruction * reconstruction,
          const Options * opts,
          float max_depth )
{
  RNTime t;
  t.Read();
  printf("Writing ply file using FETReconstruction!");
  printf("Writing file to %s\n", filename );
  int n_points = 0;
  for( int i = 0 ; i < opts->frame_end; ++i )
  {
    FETShape * shape = reconstruction->Shape(i);
    for( int j = 0 ; j < shape->NFeatures() ; ++j )
    {
      FETFeature * f = shape->Feature( j );
      R3Point p = f->Position(FALSE);
      if ( -p.Z() < max_depth ) n_points++;
    }
  }

  FILE * fid = fopen( filename, "w" );
  if ( fid )
  {
    fprintf( fid, "ply\n" );
    fprintf( fid, "format ascii 1.0\n" );
    fprintf( fid, "element vertex %d\n", n_points );
    fprintf( fid, "property float x\n" );
    fprintf( fid, "property float y\n" );
    fprintf( fid, "property float z\n" );
    fprintf( fid, "property uchar red\n" );
    fprintf( fid, "property uchar green\n" );
    fprintf( fid, "property uchar blue\n" );
    fprintf( fid, "end_header\n" );

    for( int i = 0 ; i < opts->frame_end; ++i )
    {
      FETShape * shape = reconstruction->Shape(i);
      for( int j = 0 ; j < shape->NFeatures() ; ++j )
      {
        FETFeature * f = shape->Feature( j );
        R3Point p = f->Position( FALSE );
        if ( -p.Z() < max_depth ) 
        {
          R3Point pos = f->Position( TRUE );
          RNRgb col = f->Color();
          fprintf( fid, "%f %f %f %d %d %d\n", 
                      pos.X(), pos.Y(), pos.Z(), 
                      (int)(col.R() * 255), (int)(col.G() * 255), (int)(col.B() *255) );
        }
      }
    }
  fclose( fid );
  }
  printf("Finished writing in %f\n", t.Elapsed() );
  return 0;
}

i32 io::
write_ply( const char * filename, 
          const SequenceData * seq_data,
          const Options * opts,
          float max_depth )
{
#ifdef FETREGISTER_USE_WINDOW
  FILE * f = fopen( filename, "w" );
  if ( f )
  {
    fprintf( f, "ply\n" );
    fprintf( f, "format ascii 1.0\n" );
    fprintf( f, "element vertex\n" );
    fprintf( f, "property float x\n" );
    fprintf( f, "property float y\n" );
    fprintf( f, "property float z\n" );
    fprintf( f, "property float nx\n" );
    fprintf( f, "property float ny\n" );
    fprintf( f, "property float nz\n" );
    fprintf( f, "property uchar red\n" );
    fprintf( f, "property uchar green\n" );
    fprintf( f, "property uchar blue\n" );
    fprintf( f, "end_header\n" );

    // Prepare storage
    bsc::img color_img;
    bsc::img_u16 depth_img;
    bsc::geometry_data pointcloud;
    pointcloud.n_vertices = 640 * 480;
    pointcloud.positions = (r32 *)calloc(3 * pointcloud.n_vertices, sizeof(r32));
    pointcloud.normals = (r32 *)calloc(3 * pointcloud.n_vertices, sizeof(r32));
    pointcloud.colors_a = (u8 *)calloc(3 * pointcloud.n_vertices, sizeof(u8));
    i32 n_points     = 0;

    for (i32 i = 0; i < opts->frame_end; i += 10 )
    {
      // Report status
      if ( opts->print_verbose )
      {
        if (i > 0) { printf("\r"); }
        printf("Frame %d/%d", i + 1, opts->frame_end);
        fflush(stdout);
      }

      // Prepare names
      char color_filename[512], depth_filename[512];
      sprintf(color_filename, "%s/%s", seq_data->color_dir,
              seq_data->color_names[i]);
      sprintf(depth_filename, "%s/%s", seq_data->depth_dir,
              seq_data->depth_names[i]);

      // Read in images
      color_img.read( color_filename );
      depth_img.read( depth_filename );

      // Decide if shifting
      bool bitshift = !strcmp(seq_data->dataset, "sun3d");
      
      prv_init_pointcloud( &pointcloud,
                           &color_img, 
                           &depth_img, 
                           seq_data->current_xforms[i],
                           seq_data->depth_intrinsics,
                           opts->downsample_factor, 
                           bitshift, 
                           0.01f, 
                           max_depth );

      n_points += pointcloud.n_vertices;
 
      for ( i32 j = 0 ; j < pointcloud.n_vertices ; ++j )
      {
        fprintf( f, "%f %f %f %f %f %f %d %d %d\n", 
                  pointcloud.positions[3 * j + 0],
                  pointcloud.positions[3 * j + 1],
                  pointcloud.positions[3 * j + 2],
                  pointcloud.normals[3 * j + 0],
                  pointcloud.normals[3 * j + 1],
                  pointcloud.normals[3 * j + 2],
                  pointcloud.colors_a[3 * j + 0],
                  pointcloud.colors_a[3 * j + 1],
                  pointcloud.colors_a[3 * j + 2] );
      }
    }

    fclose( f );

    // reopen the file and modify third line
    // NOTE: THIS IS BAD!
    FILE * f2;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    size_t n_lines = 0;
    f = fopen("test.ply", "r");
    f2 = fopen("test2.ply", "w");
    if (f2 == NULL)
        exit(EXIT_FAILURE);
    
    while ((read = getline(&line, &len, f)) != -1) {
        if ( n_lines != 2 ) fprintf( f2, "%s", line);
        else fprintf( f2, "element vertex %d\n", n_points );
        n_lines++;
    }

    fclose(f);
    fclose(f2);

    if (line)
        free(line);

    return 1;
  }
  printf( "Could not write to file %s\n", filename );
#endif
  return 0;
}

#endif //IO_IMPLEMENTATION
#endif //IO_H
