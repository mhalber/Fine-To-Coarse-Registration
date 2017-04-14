// FineToCoarseRegistration - fetbenchmark
//
// This is minimalistic program that is designed to read-conf files and
// calculating error based on them.

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
// Global consts
////////////////////////////////////////////////////////////////////////////////

// TODO(maciej): Possibly get rid of those and replace with std::vector?
#define MAX_STR_LENGTH        512
#define MAX_FRAME_NAME_LENGTH 64
#define MAX_ITER              32
#define MAX_N_GT_CORRS        2048

#define BSC_IMPLEMENTATION
#include "basics.h"

////////////////////////////////////////////////////////////////////////////////
// Types
////////////////////////////////////////////////////////////////////////////////
struct Correspondence
{
  char *name_A, *name_B;
  i32 idx_A, idx_B;
  bsc::vec3d point_A, point_B;
};

struct Options
{
  char * conf_filename;
  bool print_verbose;
};

struct SequenceData
{
  const char * format;

  char * dataset = NULL;
  i32 n_frames = 0;

  bsc::vec2i color_res;
  bsc::vec2i depth_res;

  bsc::mat3 depth_intrinsics;
  bsc::mat3 color_intrinsics;
  bsc::mat4 depth_extrinsics;
  bsc::mat4 color_extrinsics;

  char * depth_dir = NULL;
  char * color_dir = NULL;

  char * depth_intrinsics_filename = NULL;
  char * color_intrinsics_filename = NULL;
  char * depth_extrinsics_filename = NULL;
  char * color_extrinsics_filename = NULL;
  
  char (*depth_names)[MAX_FRAME_NAME_LENGTH] = NULL;
  char (*color_names)[MAX_FRAME_NAME_LENGTH] = NULL;

  // transformations
  bsc::mat4d * initial_xforms               = NULL;
  bsc::mat4d * current_xforms               = NULL;

  bsc::mat4d * xforms_history               = NULL;
  i32 n_xforms                              = 0;

  // ground truth correspondences
  char * correspondences_filename = NULL;
  Correspondence * gt_corrs       = NULL;
  i32 n_gt_corrs                  = 0;

  // matches
  char * pairwise_matches_filename     = NULL;
  bsc::mat4d * pairwise_xforms         = NULL;
  int * pairwise_xforms_inliers        = NULL;

  // errors
  r32 errors[ MAX_ITER ];
  r32 per_correspondence_error[MAX_N_GT_CORRS];

  // parametrization
  std::vector<double> parametrization;
};

#define IO_IMPLEMENTATION
#include "io.h"

////////////////////////////////////////////////////////////////////////////////
// Globals
////////////////////////////////////////////////////////////////////////////////
static SequenceData seq_data;
static Options opts;


////////////////////////////////////////////////////////////////////////////////
// Helpers 
////////////////////////////////////////////////////////////////////////////////
int std_comparef (const void * a, const void * b)
{
  return  (*(r32*)a > *(r32*)b) - (*(r32*)a < *(r32*)b);
}

////////////////////////////////////////////////////////////////////////////////
// Evaluation functions
////////////////////////////////////////////////////////////////////////////////

namespace eval
{

bool is_null( const bsc::mat4d & mat )
{
  r32 sum = 0.0f;
  for ( int i = 0; i < 16 ; ++i )
  {
    sum += fabs( mat.data[i] );
  }
  return sum < 1e-9;
}

void invalidate_corrs( SequenceData * seq_data,
                       const bsc::mat4d * xforms )
{
  if ( xforms == NULL ) return;

  std::vector<i32> valid_corr_ind;
  for ( i32 i = 0 ; i < seq_data->n_gt_corrs ; i++ )
  {
    const Correspondence * corr = &( seq_data->gt_corrs[i] );
    bsc::mat4d matA = xforms[ corr->idx_A ];
    bsc::mat4d matB = xforms[ corr->idx_B ];
    
    if ( is_null(matA) || is_null(matB) )
    {
      continue;
    }
    
    valid_corr_ind.push_back( i );
  }

  if ( valid_corr_ind.size() != seq_data->n_gt_corrs )
  {
    printf("Some correspondences invalid | %lu %d\n", valid_corr_ind.size(), seq_data->n_gt_corrs );
    Correspondence * new_corrs = (Correspondence*)malloc( valid_corr_ind.size() * sizeof(Correspondence) );
    for ( i32 i = 0 ; i < valid_corr_ind.size() ; ++i )
    {
      Correspondence old_corr = seq_data->gt_corrs[ valid_corr_ind[i] ];
      Correspondence new_corr;
      new_corr.name_A = strdup( old_corr.name_A );
      new_corr.name_B = strdup( old_corr.name_B );
      new_corr.idx_A = old_corr.idx_A;
      new_corr.idx_B = old_corr.idx_B;
      new_corr.point_A = old_corr.point_A;
      new_corr.point_B = old_corr.point_B;
      new_corrs[i] = new_corr;
    }

    free( seq_data->gt_corrs );
    seq_data->gt_corrs = new_corrs;
    seq_data->n_gt_corrs = valid_corr_ind.size();
    
  }
}

r32
calculate_rmse( const Correspondence * corrs,
                const i32 n_corrs,
                const bsc::mat4d * xforms,
                r32 * per_correspondence_error = NULL )
{
  r32 error = 0.0f;

  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        ++corr_idx )
  {
    const Correspondence * corr = &(corrs[corr_idx]);

    bsc::vec4d pt_A = bsc::vec4d( corr->point_A, 1.0 );
    bsc::vec4d pt_B = bsc::vec4d( corr->point_B, 1.0 );

    pt_A = xforms[corr->idx_A] * pt_A;
    pt_B = xforms[corr->idx_B] * pt_B;
    if ( is_null(xforms[corr->idx_A] ) || is_null(xforms[corr->idx_B] ) )
    {
      printf("calculate_rmse : This should never print!\n");
    }
    r32 cur_err = bsc::norm( bsc::vec3d(pt_A) - bsc::vec3d(pt_B) );

    if ( per_correspondence_error ) 
      per_correspondence_error[corr_idx] = cur_err;

    error += bsc::square(cur_err);
  }
  error = sqrt( error / (r32)n_corrs );
  return error;
}

r32
calculate_mean( const Correspondence * corrs,
                const i32 n_corrs,
                const bsc::mat4d * xforms )
{
  r32 mean_dist = 0.0f;

  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        ++corr_idx )
  {
    const Correspondence * corr = &(corrs[corr_idx]);

    bsc::vec4d pt_A = bsc::vec4d( corr->point_A, 1.0 );
    bsc::vec4d pt_B = bsc::vec4d( corr->point_B, 1.0 );

    pt_A = xforms[corr->idx_A] * pt_A;
    pt_B = xforms[corr->idx_B] * pt_B;

    mean_dist += bsc::norm( bsc::vec3d(pt_A) - bsc::vec3d(pt_B) );
  }
  
  return mean_dist / (r32)n_corrs;
} 

r32
calculate_stddev( const Correspondence * corrs,
                  const i32 n_corrs,
                  const r32 mean,
                  const bsc::mat4d * xforms )
{
  r32 squared_dist = 0.0f;

  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        ++corr_idx )
  {
    const Correspondence * corr = &(corrs[corr_idx]);

    bsc::vec4d pt_A = bsc::vec4d( corr->point_A, 1.0 );
    bsc::vec4d pt_B = bsc::vec4d( corr->point_B, 1.0 );

    pt_A = xforms[corr->idx_A] * pt_A;
    pt_B = xforms[corr->idx_B] * pt_B;

    r32 cur_dist = bsc::norm( bsc::vec3d(pt_A) - bsc::vec3d(pt_B) );
    squared_dist += cur_dist * cur_dist; 
  }

  return (squared_dist / (r32)n_corrs - mean * mean);
}

r32
calculate_median( const Correspondence * corrs,
                  const i32 n_corrs,
                  const bsc::mat4d * xforms )
{
  r32 distances[n_corrs];

  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        ++corr_idx )
  {
    const Correspondence * corr = &(corrs[corr_idx]);

    bsc::vec4d pt_A = bsc::vec4d( corr->point_A, 1.0 );
    bsc::vec4d pt_B = bsc::vec4d( corr->point_B, 1.0 );

    pt_A = xforms[corr->idx_A] * pt_A;
    pt_B = xforms[corr->idx_B] * pt_B;

    distances[corr_idx] = bsc::norm( bsc::vec3d(pt_A) - bsc::vec3d(pt_B) );
  }
  
  qsort( distances, n_corrs, sizeof(r32), std_comparef );

  r32 median = NAN;
  if (n_corrs%2==0)
  {
    median = distances[n_corrs/2] + distances[n_corrs/2+1] / 2.0f;
  }
  else
  {
    median = distances[n_corrs/2];
  }

  return median;
} 

r32
calculate_max( const Correspondence * corrs,
               const i32 n_corrs,
               const bsc::mat4d * xforms )
{
  r32 max_distance = -1e9;

  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        ++corr_idx )
  {
    const Correspondence * corr = &(corrs[corr_idx]);

    bsc::vec4d pt_A = bsc::vec4d( corr->point_A, 1.0 );
    bsc::vec4d pt_B = bsc::vec4d( corr->point_B, 1.0 );

    pt_A = xforms[corr->idx_A] * pt_A;
    pt_B = xforms[corr->idx_B] * pt_B;

    r32 distance= bsc::norm( bsc::vec3d(pt_A) - bsc::vec3d(pt_B) );
    if(distance > max_distance ) max_distance = distance;
  }
  
  return max_distance;
} 

void
calculate_error_stats( const SequenceData *seq_data,
                       float *per_correspondence_errors )
{
  float rmse = calculate_rmse( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                               seq_data->current_xforms, 
                               per_correspondence_errors );
  float mean = calculate_mean( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                               seq_data->current_xforms );
  float dev = calculate_stddev( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                                mean, seq_data->current_xforms );
  printf(" RMSE: %f, Mean %f, StdDev, %f\n", rmse, mean, dev );
}
}

void parse_arguments( int argc, char **argv, 
                      Options *opts )
{
  bsc::arg_parse args;

  args.add(bsc::argument<char *>("conf_filename",
                                 "Input configuration filename",
                                 &(opts->conf_filename)));

  args.add(bsc::argument<bool>("-v",
                               "--verbose",
                               "Prints additional information",
                               &(opts->print_verbose), 0));

  args.parse(argc, argv);
}

int main( int argc, char **argv )
{
  parse_arguments( argc, argv, &opts );
  if( !io::read_configuration_file( opts.conf_filename,
                                    &seq_data,
                                    &opts ) )
  {
    return 0;
  }
  else
  {
    float *per_correspondence_errors = NULL;
    if( opts.print_verbose ) 
    {
      per_correspondence_errors = (float*)malloc( seq_data.n_gt_corrs * sizeof(float) );
    }
    eval::calculate_error_stats( &seq_data, per_correspondence_errors );
    if( opts.print_verbose ) 
    {
      printf("  Per correspondence errors:\n");
      {
        for( int i = 0 ; i < seq_data.n_gt_corrs ; ++i )
        {
          printf("   | %4d - %6.5f\n", i, per_correspondence_errors[i] );
        }
      }
    }
  }
}