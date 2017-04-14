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
// Simple conf reading
////////////////////////////////////////////////////////////////////////////////

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

    if ( strcmp( cmd, "n_images" ) == 0 || strcmp( cmd, "n_frames" ) == 0   )
    {
      success &= parse_n_frames_cmd(buffer, cmd, seq_data );
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
    printf("\n");
  }

  // Return success
  return 1;
}

#endif //IO_IMPLEMENTATION
#endif //IO_H
