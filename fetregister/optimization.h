// FineToCoarseRegistration - fetregister - optimization
//
// This file contains all routines for gathering constraints using fine-to-coarse
// strategy, as well as optimization of that set of constraints. Functions are 
// grouped into namespaces:
//  optproxy  - planar proxy-related helper functions
//  optcorrs  - functions regarding closest point correspondences
//  optstruct - functions regarding the structural model (hierarchical and 
//              geometrical constraints)
//  opteq     - funcitons for setting up system of equations

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



#ifndef OPTIMIZATION_H_
#define OPTIMIZATION_H_

////////////////////////////////////////////////////////////////////////////////
// interface
////////////////////////////////////////////////////////////////////////////////

struct MatchCandidate
{
  int idxA;
  int idxB;
  float distance_threshold;
  float angle_threshold;
  float distance;
};

struct ThreadSafeCorrespondence
{
  FETFeature *feature1;
  FETFeature *feature2;
  double affinity;
  ThreadSafeCorrespondence( FETFeature *f1, FETFeature *f2, double a ) : feature1(f1),
                                                                         feature2(f2),
                                                                         affinity(a) {}
  ~ThreadSafeCorrespondence() {}
  ThreadSafeCorrespondence( const ThreadSafeCorrespondence &other ) { feature1 = other.feature1;
                                                                      feature2 = other.feature2;
                                                                      affinity = other.affinity; }
};


////////////////////////////////////////////////////////////////////////////////
// Proxy Helpers
////////////////////////////////////////////////////////////////////////////////
namespace optproxy
{
void get_shape_ind_connected_to_proxy( const Proxy * proxy, 
                                       std::vector<int> & ind );

std::vector<Proxy*> get_proxies_at_level( const SequenceData * seq_data, const int level );

void get_proxy_inliers( const Proxy * proxy, std::vector<FETFeature*> &inliers );

i32 get_pair_id( i32 a, i32 b, i32 size );

i32 get_pair_id_order_dependend( i32 a, i32 b, i32 size );

void reset_seen_index_increments( std::unordered_map<int, SeenInfo> * seen_count );
}


////////////////////////////////////////////////////////////////////////////////
// Closest Point Correspondences
////////////////////////////////////////////////////////////////////////////////
namespace optcorrs
{
void gather_neighbors( const int active_index, 
                       const int shape_start_idx, const int shape_end_idx,
                       std::vector<int> &neighbor_indices, int &n_neighbors,
                       float &prev_dist, float &next_dist,
                       const OptimizationParameters * optim_params,
                       const Options * opts,
                       const std::vector<double> &parametrization );

void subsample_correspondences( SequenceData * seq_data,
                                OptimizationParameters * optim_params, 
                                int n_frames );

void clear_correspondences( std::vector<FETCorrespondence*> & corrs );

double calculate_threshold( int idx_A, int idx_B, int iter_idx,
                            double max_dist, const std::vector<double> &parametrization,
                            double min, double max );

RNScalar calculate_affinity( const FETReconstruction * reconstruction,
                             const FETFeature * feature1,
                             const FETFeature * feature2 );

void create_closest_point_correspondences( FETReconstruction *reconstruction,
                                      std::vector<ThreadSafeCorrespondence> *cp_corrs,
                                      FETShape *shape1, FETShape *shape2,
                                      size_t expected_n_corrs,
                                      float max_depth,
                                      const int thread_idx,
                                      int type = COINCIDENT_RELATIONSHIP,
                                      const double *proportions = NULL );

void select_closest_point_correspondences( const std::vector<ThreadSafeCorrespondence> *input_corrs,
                                      std::vector<ThreadSafeCorrespondence> *output_corrs,
                                      int & output_idx, int min_idx, int max_idx,
                                      int n_corrs, const int thread_idx, const int active_idx );

void add_closest_point_correspondeces( FETReconstruction *reconstruction,
                                       SequenceData *seq_data,
                                       OptimizationParameters *optim_params,
                                       const Options *opts,
                                       const bool modify_weight = true );
}


////////////////////////////////////////////////////////////////////////////////
// Structural and To Model Correspondences
////////////////////////////////////////////////////////////////////////////////

namespace optstruct
{

void add_relationship_between_proxies( Proxy * proxy_A,
                                  Proxy * proxy_B,
                                  SequenceData * seq_data,
                                  const OptimizationParameters * optim_params );

void create_geometrical_relationships( int level_idx, 
                                   SequenceData *seq_data,
                                   const OptimizationParameters *optim_params );

void create_hierarchical_structure( int level_idx, int end_frame,
                                    SequenceData * seq_data,
                                    OptimizationParameters * optim_params,
                                    const Options * opts );

void add_geometrical_correspondeces( const std::vector<StructRel> *relationships, 
                                     std::vector<FETCorrespondence*> & corrs,
                                     const Options * opts );


void add_hierarchical_correspondences( const FETReconstruction * reconstruction,
                                       const StructuralModel * model,
                                       std::vector<FETCorrespondence*> & corrs,
                                       const Options * opts );

}


////////////////////////////////////////////////////////////////////////////////
// Equations
////////////////////////////////////////////////////////////////////////////////
namespace opteq
{

void add_trajectory_equations( FETReconstruction *reconstruction,
                               const SequenceData *data,
                               RNSystemOfEquations *system,
                               const OptimizationParameters *optim_params,
                               const Options *opts );

void add_correspondence_equations( FETReconstruction *reconstruction,
                                   RNSystemOfEquations *system,
                                   const std::vector<FETCorrespondence*> &corrs,
                                   RNScalar total_weight,
                                   const RNScalar * input_proportions,
                                   std::string type_string,
                                   const Options *opts );

}

////////////////////////////////////////////////////////////////////////////////
// Main Optimization
////////////////////////////////////////////////////////////////////////////////

void optimize_transformations( FETReconstruction *reconstruction,
                               const SequenceData *seq_data,
                               const OptimizationParameters *optim_params,
                               const Options *opts);

#endif

#ifdef OPTIMIZATION_IMPLEMENTATION


////////////////////////////////////////////////////////////////////////////////
// Proxy Helpers
////////////////////////////////////////////////////////////////////////////////
namespace optproxy
{
void
get_shape_ind_connected_to_proxy( const Proxy * proxy, std::vector<int> & ind )
{
  if( proxy->level == 0 )
  {
    if( !proxy->parent && proxy->children.NEntries() == 0 ) return;

    FETFeature * feature = proxy->feature;
    FETCorrespondence * corr = feature->Correspondence( 0 );
    if( !corr ) return;
    
    FETFeature * f1 = corr->Feature( 0 );
    FETFeature * f2 = corr->Feature( 1 );
    if( f1 && f2 ) // valid corr
    {
      FETShape * shape = f1->Shape();
      if( shape ) ind.push_back( shape->reconstruction_index );
    }
    return;
  }
  else
  {
    for ( int child_idx = 0 ; child_idx < proxy->children.NEntries() ; ++child_idx )
    {
      Proxy * child = proxy->children[ child_idx ];
      get_shape_ind_connected_to_proxy( child, ind );
    }
  }
}

std::vector<Proxy*>
get_proxies_at_level( const SequenceData * seq_data, const int level )
{
  Structure * structure = seq_data->model->structure;
  std::vector<Proxy*> proxies_at_level;
  for( int proxy_idx = 0 ; proxy_idx < structure->NProxies() ; ++proxy_idx )
  {
    Proxy * proxy = structure->GetProxy( proxy_idx );
    if( proxy && proxy->level == level ) proxies_at_level.push_back( proxy );
  }
  return proxies_at_level;
}

void
get_proxy_inliers( const Proxy * proxy, std::vector<FETFeature*> &inliers )
{

  if( proxy->level > 0 )
  {
    for( int i = 0 ; i < proxy->children.NEntries() ; ++i )
    {
      get_proxy_inliers( proxy->children.Kth( i ), inliers );
    }
  } 
  else
  {
    for( int i = 0 ; i < proxy->shape->NFeatures() ; ++i )
    {
      FETFeature * shape_feature = proxy->shape->Feature(i);
      for( int j = 0 ; j < shape_feature->NCorrespondences() ; ++j )
      {
        FETCorrespondence * corr = shape_feature->Correspondence(j);
        inliers.push_back( corr->Feature(1) );
      }
    }
  }
}

i32 
get_pair_id( i32 a, i32 b, i32 size )
{
  i32 arr[2] = {a, b};
  qsort(arr, 2, sizeof( i32 ), std_compare );
  return arr[1] * size + arr[0];
}

i32 
get_pair_id_order_dependend( i32 a, i32 b, i32 size )
{
  return a * size + b;
}

void
reset_seen_index_increments( std::unordered_map<int, SeenInfo> * seen_count )
{
  std::unordered_map<int,SeenInfo>::iterator seen_it;
  for ( seen_it = seen_count->begin();
        seen_it != seen_count->end(); 
        seen_it++ )
  {
    seen_it->second.incremented = false;
  }
}
}

////////////////////////////////////////////////////////////////////////////////
// Correspondence Helpers 
////////////////////////////////////////////////////////////////////////////////
namespace optcorrs
{
void
gather_neighbors( const int active_index, 
                  const int shape_start_idx, const int shape_end_idx,
                  std::vector<int> &neighbor_indices, int &n_neighbors,
                  float &prev_dist, float &next_dist,
                  const OptimizationParameters * optim_params,
                  const Options * opts,
                  const std::vector<double> &parametrization )
{
  float active_dist = parametrization[active_index];
  // walk back 
  if( active_index > shape_start_idx )
  {
    int prev_idx = active_index - 1;
    prev_dist = parametrization[ prev_idx ];
    while( fabs( active_dist-prev_dist ) < 0.5*optim_params->segment_length )
    {
      neighbor_indices[ n_neighbors ] = prev_idx;
      n_neighbors++; prev_idx--;
      
      if( prev_idx >= 0 ) prev_dist = parametrization[ prev_idx ];
      else                break;

      if( n_neighbors >= opts->frame_end ) break;
    }
  }

  // walk forward
  if( active_index < shape_end_idx-1 )
  {
    int next_idx = active_index + 1;
    next_dist = parametrization[ next_idx ];
    while( fabs( active_dist-next_dist ) < 0.5*optim_params->segment_length )
    {
      neighbor_indices[ n_neighbors ] = next_idx;
      n_neighbors++; next_idx++;

      if( next_idx <= opts->frame_end - 1 ) next_dist = parametrization[ next_idx ];
      else                                    break;
    
      if( n_neighbors >= opts->frame_end ) break;
    }
  }
}

void 
subsample_correspondences( SequenceData * seq_data,
                           OptimizationParameters * optim_params, 
                           int n_frames )
{
  FETReconstruction * reconstruction = seq_data->reconstruction;
  i32 n_corrs = optim_params->max_n_corrs * 0.5f;
  if( reconstruction->NCorrespondences() <= n_corrs ) return; 

  long long n_tries = 0;
  while (reconstruction->NCorrespondences() > n_corrs) 
  {
    FETCorrespondence * corr = reconstruction->Correspondence( get_random_int(0, 0, reconstruction->NCorrespondences()-1));
    FETFeature * f = corr->Feature(1);
    if( f->NCorrespondences() >= 5 )
    {
      delete corr;
    }

    n_tries++;
    if(n_tries > 1e12 ) break;
  }
}

void
clear_correspondences( std::vector<FETCorrespondence*> & corrs )
{
  for ( i32 corr_idx = 0 ; 
        corr_idx < corrs.size() ;
        corr_idx++ )
  {
    delete corrs[corr_idx];
  }
  corrs.clear();
}

double calculate_threshold( int idx_A, int idx_B, int iter_idx,
                            double max_dist, const std::vector<double> &parametrization,
                            double min, double max )
{
  // base function
  double e = 0.5;
  double a = 1.0 / std::pow(max_dist, e);
  double d = std::fabs( parametrization[idx_A] - parametrization[idx_B] );
  double t = a * std::pow(d, e);

  return std::min( max, min + t * (max - min) );
}

}

namespace optcorrs
{
////////////////////////////////////////////////////////////////////////////////
// Closest Point Correspondences
////////////////////////////////////////////////////////////////////////////////

double
calculate_affinity( const FETReconstruction * reconstruction,
                  const FETFeature * feature1,
                  const FETFeature * feature2 )
{
  // Get/check shapes
  if( !feature1 || !feature2 ) return -1.0f;
  FETShape *shape1 = feature1->shape;
  FETShape *shape2 = feature2->shape;
  if( !shape1 || !shape2) return -1.0f;
  if( feature1->GeneratorType() != feature2->GeneratorType() ) return -1.0f;
  if( feature1->ShapeType() != feature2->ShapeType() ) return -1.0f;

  // Keep going back and forth on this
  double affinity = fabs(0.5 * (feature1->Salience() + feature2->Salience()));

  R3Vector n1 = feature1->Normal( TRUE );
  R3Vector n2 = feature2->Normal( TRUE );
  n1.Normalize();
  n2.Normalize();
  double dot = n1.Dot(n2);
  double normal_angle = (dot < 1) ? acos(dot) : 0;

  if( feature1->ShapeType() == PLANE_FEATURE_SHAPE )
  {
    if( normal_angle > reconstruction->max_normal_angle ) return -1.0f; // normals deviate too much.
    affinity *= (1.0 - 0.36 * (normal_angle / reconstruction->max_normal_angle));
  }
  
  else if( feature1->ShapeType() == LINE_FEATURE_SHAPE && feature2->ShapeType() == LINE_FEATURE_SHAPE )
  { 
    R3Vector d1 = feature1->Direction( TRUE );
    R3Vector d2 = feature2->Direction( TRUE );
    d1.Normalize();
    d2.Normalize();
    double dot = fabs(d1.Dot(d2));
    double direction_angle = acos( std::min(dot, 1.0) );
    affinity *= (1.0 - 0.25 * (direction_angle / reconstruction->max_normal_angle));

    if( direction_angle > reconstruction->max_normal_angle ) return -1.0f;

    if( (feature1->GeneratorType() == VALLEY_FEATURE_TYPE && feature2->GeneratorType() == VALLEY_FEATURE_TYPE) ||
         (feature1->GeneratorType() == RIDGE_FEATURE_TYPE && feature2->GeneratorType() == RIDGE_FEATURE_TYPE)  )
    {
      if( normal_angle > reconstruction->max_normal_angle ) return -1.0f;
      else
      {
        affinity *= (1.0 - 0.25 * (normal_angle / reconstruction->max_normal_angle));
      }
    }
    else if( feature1->GeneratorType() == SILHOUETTE_FEATURE_TYPE )
    {
      if( normal_angle > (0.5 * RN_PI) ) return -1.0f;
      else
      {
        affinity *= (1.0 - 0.25 * (normal_angle / (0.5 * RN_PI)));
      }
    } 
    else
    {
      return -1.0f;
    }
  }
  else
  {
    return -1.0f;
  }

  return affinity;
}

void 
create_closest_point_correspondences( FETReconstruction *reconstruction,
                                      std::vector<ThreadSafeCorrespondence> *cp_corrs,
                                      FETShape *shape1, FETShape *shape2,
                                      size_t expected_n_corrs,
                                      float max_depth,
                                      const int thread_idx,
                                      int type,
                                      const double *proportions )
{
  // Check number of correspondences
  if(expected_n_corrs <= 0)
  { 
    return;
  }

  // Variables
  double total_salience = 0;
  std::vector<FETFeature*> query_features; 
  std::vector<bool> selected( shape2->NFeatures(), false ); 

  size_t n_valid_features = 0;
  for ( int i = 0 ; i < shape2->NFeatures() ; ++i )
  {
    FETFeature *feature = shape2->Feature( i );
    if( proportions && 
         proportions[feature->GeneratorType()] == 0.0 ) continue;
    if( (feature->Salience() == 0.0f ) || 
         (feature->IsOnBoundary() &&
         (!feature->IsOnSilhouetteBoundary()))) continue;
    n_valid_features++;
  }
  // if( shape1->reconstruction_index == 3150 && shape2->reconstruction_index == 3090 ) 
    // printf("          -- STATS : %d %f \n", n_valid_features, expected_n_corrs ); 

  // select features
  query_features.reserve( expected_n_corrs );
  if( n_valid_features < expected_n_corrs )
  {  
    // if( shape1->reconstruction_index == 3150 && shape2->reconstruction_index == 3090 ) 
    // printf("          -- JUST COPYING\n"); 
    for ( int i = 0 ; i < shape2->NFeatures() ; ++i )
    {
      FETFeature *feature = shape2->Feature( i );
      if( proportions && 
           proportions[feature->GeneratorType()] == 0.0 ) {continue;}
      if( (feature->Salience() == 0.0f ) || 
           (feature->IsOnBoundary() &&
           (!feature->IsOnSilhouetteBoundary()))) {continue;}

      query_features.push_back( feature );
      total_salience += feature->Salience();
    }
  }
  else {
    //  if( shape1->reconstruction_index == 3150 && shape2->reconstruction_index == 3090 ) 
      // printf("          -- SAMPLING! %lu\n", expected_n_corrs ); 
    // select query features randomly
    while( query_features.size() < expected_n_corrs )
    {
      int feat_idx = get_random_int( thread_idx, 0, shape2->NFeatures() - 1 );
      FETFeature *feature2 = shape2->Feature( feat_idx );
      
      // validity tests
      if( !feature2 ) continue;
      if( selected[ feat_idx ] ) continue;
      if( proportions && 
          proportions[feature2->GeneratorType()] == 0.0 ) continue;
      if( (feature2->Salience() == 0.0f ) || 
          (feature2->IsOnBoundary() &&
          (!feature2->IsOnSilhouetteBoundary()))) continue;

      // random exclusion
      double salience = feature2->Salience();
      if( get_random_real( thread_idx ) > salience ) continue;

      query_features.push_back( feature2 );
      selected[ feat_idx ] = true;
      total_salience += salience; 
    }
    //  printf("          -- DONE SAMPLING! %lu\n", expected_n_corrs ); 
  }

  // Compute correspondences for shape2 -> shape1
  for (int i = 0; i < query_features.size(); i++)
  {
    FETFeature *feature2 = query_features[i];
    if( !feature2 ) continue;


    // Find closest feature on shape1
    FETFeature *feature1 = shape1->FindClosestFeature(
        feature2, shape2->current_transformation,
        0.0,
        reconstruction->max_euclidean_distance,
        reconstruction->max_descriptor_distances,
        reconstruction->max_normal_angle,
        reconstruction->min_distinction,
        reconstruction->min_salience,
        reconstruction->discard_boundaries ); 

    if( !feature1 )
    {
      continue;
    }
    
    if( feature1->Position().Z() < -max_depth || 
         feature2->Position().Z() < -max_depth )
    {
      continue;
    }

    // check ifhave same shape
    if( feature1->ShapeType() != feature2->ShapeType() )
    {
      continue;
    }

    // check ifhave same type..
    if( feature1->GeneratorType() != feature2->GeneratorType() )
    {
      continue;
    }

    // Check ifshould discard because on boundary
    if( reconstruction->discard_boundaries && feature1->IsOnBoundary() )
    {
      if(!feature1->IsOnSilhouetteBoundary() ||
          !feature2->IsOnSilhouetteBoundary())
      {
        continue;
      }
    }

    // Compute/check affinity
    RNScalar affinity = calculate_affinity( reconstruction, feature1, feature2 );

    // Create correspondence with sorted shapes
    if( shape1->reconstruction_index < shape2->reconstruction_index )
    {
      if( affinity > 0.0 ) cp_corrs->push_back( ThreadSafeCorrespondence( feature1, feature2, affinity ) );
    }
    else
    {
      if( affinity > 0.0) cp_corrs->push_back( ThreadSafeCorrespondence( feature2, feature1, affinity ) );
    }
  }
  // if( shape1->reconstruction_index == 3150 && shape2->reconstruction_index == 3090 ) 
      // printf("          -- ALL DONE! \n" ); 
}

void 
select_closest_point_correspondences( const std::vector<ThreadSafeCorrespondence> *input_corrs,
                                   std::vector<ThreadSafeCorrespondence> *output_corrs,
                                   int & output_idx, int min_idx, int max_idx,
                                   int n_corrs, const int thread_idx, const int active_idx )
{
 
  /* Check trivial conditions */
  if( n_corrs < 0)
    return;
  if( input_corrs->empty() )
    return;
  
// #ifdef VISUALIZE_STATS
//   save_per_shape_stats( input_corrs, n_shapes, grid_a );
// #endif

  std::vector<ThreadSafeCorrespondence> valid_corrs;
  valid_corrs.reserve( input_corrs->size() );
  for(int i = 0 ; i < input_corrs->size() ; ++i )
  {
    if( (*input_corrs)[i].feature1 != NULL && (*input_corrs)[i].feature2 != NULL )
    {
      valid_corrs.push_back( (*input_corrs)[i] );
    }
  }
  /* ifnot too many, then copy all */
  if( valid_corrs.size() < n_corrs )
  {
    for( int i = 0 ; i < valid_corrs.size() ; ++i )
    {
      if( (output_idx >= max_idx) || (output_idx < min_idx) ) break;
      (*output_corrs)[output_idx] = valid_corrs[i]; 
      output_idx++;
    }
    return;
  }
  /* otherwise we need to select some */
  else
  {
    int n_valid_corrs = valid_corrs.size();
    int n_selected_corrs = 0;

    std::vector<bool> select( valid_corrs.size(), false );
    // int n_tries = 0;
    while( n_selected_corrs < n_corrs )
    {
      if( (output_idx >= max_idx) || (output_idx < min_idx) ) break;
      int idx = get_random_int( thread_idx, 0, n_valid_corrs - 1 );

      if( !select[idx] )
      {
        (*output_corrs)[output_idx] = valid_corrs[idx];
        select[idx] = true;
        output_idx++;
        n_selected_corrs++;
      }
      // n_tries++;
    }
  }

// #ifdef VISUALIZE_STATS
//   save_per_shape_stats( output_corrs, n_shapes, grid_b );
// #endif

}

int
InternalDiscardOutliers(std::vector<ThreadSafeCorrespondence> &cp_corrs,
                        RNScalar max_zscore)
{
  // Check correspondences
  if(max_zscore >= 10)
    return cp_corrs.size();
  if(cp_corrs.size() < 5)
    return cp_corrs.size();

  // Allocate distances
  RNScalar *distances = new RNScalar[cp_corrs.size()];

  // Sum distances
  RNScalar sum = 0;
  RNScalar n_valid = 0;
  for (int i = 0; i < cp_corrs.size(); i++)
  {
    ThreadSafeCorrespondence &correspondence = cp_corrs[i];
    if(correspondence.feature1 && correspondence.feature2 )
    {
      distances[i] = (correspondence.feature1->Position(TRUE) - correspondence.feature2->Position(TRUE)).Length();
      sum += distances[i];
      n_valid++;
    }
  }

  // Compute mean
  RNScalar mean = sum / n_valid;

  // Sum squared residuals
  RNScalar rss = 0;
  for (int i = 0; i < cp_corrs.size(); i++)
  {
    if(cp_corrs[i].feature1 && cp_corrs[i].feature2 )
    {
      RNScalar delta = distances[i] - mean;
      rss += delta * delta;
    }
  }

  // Compute minimum affinity
  RNScalar stddev = sqrt(rss / n_valid);
  RNScalar max_distance = mean + max_zscore * stddev;

  // Mark outliers
  int n_outliers = 0;
  for (int i = 0; i < cp_corrs.size(); i++)
  {
    if( distances[i] > max_distance )
    {
      ThreadSafeCorrespondence &corr = cp_corrs[i];
      corr.feature1=NULL;
      corr.feature2=NULL;
      n_outliers++;
    }
  }

  // Delete distances
  delete[] distances;

  // return how many corrs left
  return n_valid - n_outliers;
}

void 
DiscardOutlierCorrespondences(std::vector<ThreadSafeCorrespondence> &cp_corrs,
                              RNScalar max_zscore = 3, 
                              int max_iterations = 8, 
                              int min_correspondences = 5)
{
  // Iteratively discard correspondences until convergence
  int n = cp_corrs.size();
  for (int i = 0; i < max_iterations; i++)
  {
    if(n <= min_correspondences)
      break;
    n = InternalDiscardOutliers(cp_corrs, max_zscore);
  }
}

void add_closest_point_correspondences_thread( SequenceData *seq_data,
                                               std::vector<ThreadSafeCorrespondence> *output_corrs,
                                               OptimizationParameters *optim_params,
                                               const Options *opts,
                                               const int thread_idx,
                                               const int shape_start_idx, 
                                               const int shape_end_idx,
                                               const int corr_start_idx,
                                               const int corr_end_idx,
                                               const float total_allocation,
                                               const int max_n_corrs,
                                               const std::vector<double> parametrization )
{
  /* setup helper storage */
  FETReconstruction * reconstruction = seq_data->reconstruction;
  int corr_idx = corr_start_idx;

  if( opts->print_verbose ) printf("      -- Launching thread %d\n", thread_idx );

  /* add the correspondences */
  for( int j = shape_start_idx ; j != shape_end_idx; ++j )
  {
    // printf("        -- Index : %d\n", j ); 
    int active_index = j;
    FETShape * active_shape = reconstruction->Shape( active_index );
    if( !active_shape ) continue;
    if( active_shape->NFeatures() == 0 ) continue;

    R3Box active_bbox = seq_data->bboxes[ active_index ];
    // active_bbox.Inflate( 1.25f );

    /* clear data */
    std::vector<int> neighbor_indices( opts->frame_end, 0 );
    int n_neighbors = 0;
    std::vector<ThreadSafeCorrespondence> cp_corrs_tmp;

    float prev_dist = 0.0f, next_dist = 0.0f;
    gather_neighbors( active_index, shape_start_idx, shape_end_idx,
                      neighbor_indices, n_neighbors, 
                      prev_dist, next_dist,
                      optim_params, 
                      opts,
                      parametrization );

    /* given neighbors, attempt to create correspondences */
    float cur_allocation = next_dist - prev_dist;
    for( int k = 0 ; k < n_neighbors ; ++k )
    {
      int neighbor_index = neighbor_indices[ k ];

      FETShape *neighbor_shape = reconstruction->Shape( neighbor_index );

      if( !neighbor_shape ) {continue;}
      if( neighbor_shape->NFeatures() == 0 ) {continue;}

      R3Box neighbor_bbox = seq_data->bboxes[ neighbor_index ];
      R3Box intersection = R3null_box;
      bool intersects = R3Intersects( active_bbox, neighbor_bbox, 
                                      &intersection);

      if( intersects && 
           intersection.Volume() > 0.1 * std::min( active_bbox.Volume(), neighbor_bbox.Volume() ) )
      {

        double min_dt = std::min( optim_params->init_distance_threshold, 
                                  optim_params->final_distance_threshold );
        double max_dt = std::max( optim_params->init_distance_threshold, 
                                  optim_params->final_distance_threshold );
        double dt = calculate_threshold( active_index, neighbor_index, 
                                         optim_params->iter_idx,
                                         0.5 * optim_params->segment_length, 
                                         parametrization,
                                         min_dt, max_dt );
        double min_at = std::min( optim_params->init_angle_threshold, 
                                  optim_params->final_angle_threshold );
        double max_at = std::max( optim_params->init_angle_threshold, 
                                  optim_params->final_angle_threshold );
        double at = calculate_threshold( active_index, neighbor_index,
                                         optim_params->iter_idx, 
                                         0.5 * optim_params->segment_length,
                                         parametrization,
                                         min_at, max_at );

        reconstruction->max_euclidean_distance = dt;
        reconstruction->max_normal_angle       = at;

        int n_corrs = cp_corrs_tmp.size();

        create_closest_point_correspondences( reconstruction, 
                                              &cp_corrs_tmp,
                                              neighbor_shape, active_shape, 
                                              optim_params->n_samples,
                                              optim_params->max_depth,
                                              thread_idx,
                                              COINCIDENT_RELATIONSHIP,
                                              optim_params->cp_proportions );

        if( (int)cp_corrs_tmp.size() - n_corrs > 0 )
        {
          if( dt > optim_params->glob_max_dist ) optim_params->glob_max_dist = dt;
          if( dt < optim_params->glob_min_dist ) optim_params->glob_min_dist = dt;
          if( at > optim_params->glob_max_angle ) optim_params->glob_max_angle = at;
          if( at < optim_params->glob_min_angle ) optim_params->glob_min_angle = at;
        }
      }
    }

    double percentage = ((double)cur_allocation / (double)total_allocation);
    int n_corrs = round( percentage * max_n_corrs );

    if(1)
    {
      DiscardOutlierCorrespondences(cp_corrs_tmp, 6);
    }

    // NOTE: This might be unsynchronized...
    select_closest_point_correspondences( &cp_corrs_tmp, 
                                          output_corrs, 
                                          corr_idx, 
                                          corr_start_idx, corr_end_idx, 
                                          n_corrs, thread_idx, active_index );
  }
  if( opts->print_verbose )
  {
    printf("      -- Thread %2d Finished! | Shapes %5d to %5d (%5d) | "
           "Created %6d out of %6d allocated correspondences.\n", 
           thread_idx,
           shape_start_idx, shape_end_idx, opts->frame_end,
           corr_idx - corr_start_idx,
           corr_end_idx - corr_start_idx );
    fflush(stdout);
  }
}

void
add_closest_point_correspondeces( FETReconstruction *reconstruction,
                                  SequenceData *seq_data,
                                  OptimizationParameters *optim_params,
                                  const Options *opts,
                                  const bool modify_weight )
{
  RNTime t; t.Read();
  /* preallocate storage */
  int n_threads = optim_params->n_threads;

  std::vector<int> neighbor_indices( opts->frame_end, 0 );
  int n_neighbors = 0;
  int n_shapes_per_thread = round( opts->frame_end / (float)n_threads );
  float total_allocation = 0;
  std::vector<float> allocations_per_thread(n_threads, 0);

  optim_params->glob_max_angle = -1e9;
  optim_params->glob_min_angle = 1e9;
  optim_params->glob_max_dist = -1e9;
  optim_params->glob_min_dist = 1e9;

  /* report useful info */
  if( opts->print_verbose )
  {
    printf("   Creating closest point correspondences\n"
          "      -- Segment Length:     %5.3f / %5.3f\n"
          "      -- Distance Threshold: %5.3f m\n"
          "      -- Angle Threshold:    %5.3f rad(%f deg)\n",
              optim_params->segment_length, 
              seq_data->parametrization[opts->frame_end-1],
              optim_params->distance_threshold,
              optim_params->angle_threshold,
              bsc::rad2deg( optim_params->angle_threshold ) );
  }

  /* for each shape, gather neighbors and compute correspondence allocation */
  for( int i = 0 ; i < n_threads ; ++i )
  {
    int start_idx = i * n_shapes_per_thread;
    int end_idx   = std::min( opts->frame_end, (i+1) * n_shapes_per_thread);
    if( i > 0 ) allocations_per_thread[i] = allocations_per_thread[i-1];
    for( int j = start_idx ; j != end_idx; ++j )
    {
      int active_index = j;
      FETShape * active_shape = reconstruction->Shape(active_index);
      if( !active_shape ) continue;
      if( active_shape->NFeatures() == 0 ) continue;

      R3Box active_bbox = seq_data->bboxes[ active_index ];
      active_bbox.Inflate( 1.25f );

      /* clear neighbor info */
      for( int k = 0; k < opts->frame_end; ++k) neighbor_indices[k] = 0;
      n_neighbors = 0;

      float prev_dist = 0.0f, next_dist = 0.0f;
      gather_neighbors( active_index, start_idx, end_idx,
                        neighbor_indices, n_neighbors, 
                        prev_dist, next_dist,
                        optim_params, 
                        opts,
                        seq_data->parametrization );

      total_allocation += next_dist - prev_dist;
      allocations_per_thread[i] += next_dist - prev_dist;
    }
  }

  /* according to volume allocation per thread, compute starting indices into storage */
  std::vector<int> thread_start_indices;
  thread_start_indices.push_back(0);
  for ( int i = 0 ; i < allocations_per_thread.size(); ++i )
  {
    double frac = (allocations_per_thread[i]/total_allocation);
    thread_start_indices.push_back( frac * optim_params->max_n_corrs );
  }

  /* for each shape, gather neighbors and compute correspondences */
  std::vector< std::thread > threads;
  std::vector< ThreadSafeCorrespondence > output_corrs;
  ThreadSafeCorrespondence tcs = { NULL, NULL, 0 };
  output_corrs.resize( optim_params->max_n_corrs, tcs );
  if( opts->print_verbose ) printf("    Launching threads: \n");
  for( int i = 0 ; i < n_threads ; ++i )
  {
    int shape_start_idx = i * n_shapes_per_thread;
    int shape_end_idx   = std::min( opts->frame_end, (i+1) * n_shapes_per_thread);
    int corr_start_idx  = thread_start_indices[i];
    int corr_end_idx    = thread_start_indices[i+1];


    threads.push_back( std::thread( optcorrs::add_closest_point_correspondences_thread,
                                    seq_data,
                                    &output_corrs,
                                    optim_params,
                                    opts,
                                    i,
                                    shape_start_idx, 
                                    shape_end_idx,
                                    corr_start_idx,
                                    corr_end_idx,
                                    total_allocation,
                                    optim_params->max_n_corrs,
                                    seq_data->parametrization ) );
  }

  for( int i = 0 ; i < n_threads ; i++ )
  {
    threads[i].join();
  }
  if( opts->print_verbose ) printf("    -- Threads Finished!\n");

  /* Convert ThreadSafeCorrespondences to FETCorrespondences */
  seq_data->cp_corrs.reserve( output_corrs.size() );
  long corr_type_counter[ NUM_FEATURE_TYPES ] = {};
  for ( int i = 0 ; i < output_corrs.size() ; ++i )
  {
    ThreadSafeCorrespondence &tsc = output_corrs[i];
    if( tsc.feature1 != NULL && tsc.feature2 != NULL && tsc.affinity > 1e-6 )
    {
      seq_data->cp_corrs.push_back( new FETCorrespondence( NULL, 
                                                    tsc.feature1, 
                                                    tsc.feature2, 
                                                    tsc.affinity, 
                                                    COINCIDENT_RELATIONSHIP ) );
      corr_type_counter[ tsc.feature1->GeneratorType() ]++;
    }
    else
    {
      seq_data->cp_corrs.push_back( NULL );
    }
  }
  int n_valid = std::count_if( seq_data->cp_corrs.begin(), 
                               seq_data->cp_corrs.end(), 
                               [](FETCorrespondence *c){ return c != NULL; });

  int total = std::count_if( seq_data->cp_corrs.begin(), 
                               seq_data->cp_corrs.end(), 
                               [](FETCorrespondence *c){ return c != NULL; });


  if( opts->print_verbose )
  {
    printf("    Created %d closest point correspondences in %f s.\n"
          "    Total with loop clusre corrs. is %d\n"
          "      -- Distance Thresholds:     %5.3f - %5.3f\n"
          "      --    Angle Thresholds:     %5.3f - %5.3f\n", 
                                                n_valid, t.Elapsed(), total,
                                                optim_params->glob_min_dist, 
                                                optim_params->glob_max_dist,
                                                optim_params->glob_min_angle, 
                                                optim_params->glob_max_angle );
    printf("         -- Planar      : %lu\n", 
                                         corr_type_counter[PLANE_FEATURE_TYPE]);
    printf("         -- Ridges      : %lu\n", 
                                         corr_type_counter[RIDGE_FEATURE_TYPE]);
    printf("         -- Valleys     : %lu\n", 
                                        corr_type_counter[VALLEY_FEATURE_TYPE]);
    printf("         -- Silhouettes : %lu\n", 
                                    corr_type_counter[SILHOUETTE_FEATURE_TYPE]);
  }
}
}

namespace optstruct
{
////////////////////////////////////////////////////////////////////////////////
// Structural Correspondences
////////////////////////////////////////////////////////////////////////////////

void
add_relationship_between_proxies( Proxy * proxy_A,
                                  Proxy * proxy_B,
                                  SequenceData * seq_data,
                                  const OptimizationParameters * optim_params )
{
  if( proxy_A == NULL ) { printf("   Proxy A invalid!\n"); return;}
  if( proxy_B == NULL ) { printf("   Proxy B invalid!\n"); return;}
  if( proxy_A->feature == NULL ) { printf("   Feature A invalid!\n"); return;}
  if( proxy_B->feature == NULL ) { printf("   Feature B invalid!\n"); return;}
  R3Vector n_A = proxy_A->feature->Normal( TRUE ); 
  R3Vector n_B = proxy_B->feature->Normal( TRUE );
  R3Point p_A = proxy_A->feature->Position( TRUE );
  R3Point p_B = proxy_B->feature->Position( TRUE );

  R3Plane plane_A( p_A, n_A );
  R3Plane plane_B( p_B, n_B );
  
  StructRel sr;
  sr.proxy_A = proxy_A;
  sr.proxy_B = proxy_B;

  n_A.Normalize();
  n_B.Normalize();
  double dot = n_A.Dot(n_B);
  double angle = (dot < 1) ? acos(dot) : 0;
  double angle_sigma = bsc::deg2rad(7.5); // 15 degrees of tolerance; (95% of observations 2 sigma away)
  if( fabs( RN_PI_OVER_TWO - angle ) < optim_params->structure_angle_threshold )
  {
    sr.relationship_type = PERPENDICULAR_RELATIONSHIP;
    sr.score  = mmath::gauss( angle, RN_PI_OVER_TWO, angle_sigma );
    if( sr.score > 0.1 ) seq_data->model->relationships->push_back( sr );
  }

  if( fabs( angle ) < optim_params->structure_angle_threshold )
  {
    sr.relationship_type = PARALLEL_RELATIONSHIP;
    sr.score  = mmath::gauss( angle, 0.0, angle_sigma );
    if( sr.score > 0.1 ) seq_data->model->relationships->push_back( sr );
  }

  if( fabs( RN_PI - angle ) < optim_params->structure_angle_threshold )
  {
    sr.relationship_type = ANTIPARALLEL_RELATIONSHIP;
    sr.score  = mmath::gauss( angle, RN_PI, angle_sigma );
    if( sr.score > 0.1 ) seq_data->model->relationships->push_back( sr );
  }
}

void 
create_geometrical_relationships( int level_idx, 
                               SequenceData * seq_data,
                               const OptimizationParameters * optim_params )
{
  std::vector<Proxy*> cur_proxies = optproxy::get_proxies_at_level( seq_data, level_idx );
  if( cur_proxies.empty() ) return;

  for ( int i = 0 ; i < cur_proxies.size() - 1 ; i++ )
  {
    Proxy * pA = cur_proxies[i];
    for ( int j = i + 1; j < cur_proxies.size() ; j++ )
    {
      Proxy * pB = cur_proxies[j];
      int max_lim = std::max( pA->max_shape_index, pB->max_shape_index );
      int min_lim = std::min( pA->min_shape_index, pB->min_shape_index );
      double dist = seq_data->parametrization[max_lim] -
                    seq_data->parametrization[min_lim];
      assert( dist < 0 );
      if( dist < 2.0 * optim_params->segment_length + 0.1 ) /* add a little bit of bias */
          add_relationship_between_proxies( pA, pB, seq_data, optim_params );
    }
  }
}

void
create_hierarchical_structure( int level_idx, int end_frame,
                               SequenceData * seq_data,
                               OptimizationParameters * optim_params,
                               const Options * opts )
{
  if( level_idx == 0 )
  {
    // Clear exisiting data
    if( seq_data->model )
    {
      Structure * structure = seq_data->model->structure;
      std::unordered_map<FETShape*, Proxy*> *s2p = seq_data->model->shape_to_proxy;
      std::vector<StructRel> * relationships = seq_data->model->relationships;
      
      if( structure ) { delete structure; seq_data->model->structure = NULL; }
      if( s2p ) { delete s2p ; seq_data->model->shape_to_proxy = NULL; }
      if( relationships ) { delete relationships; seq_data->model->relationships = NULL; }
      free( seq_data->model );
    }

    // Create fresh storage
    seq_data->model = (StructuralModel*)malloc( sizeof(StructuralModel));
    seq_data->model->structure      = new Structure( seq_data->reconstruction );
    seq_data->model->shape_to_proxy = new std::unordered_map<FETShape*, Proxy*>();
    seq_data->model->relationships  = new std::vector<StructRel>();
    seq_data->model->structure->CreateBaseLevelProxies( optim_params->max_depth );
    update_shape_proxy_map( seq_data->model );

    // Only use subset
    optcorrs::subsample_correspondences( seq_data, 
                                         optim_params, end_frame );
  }
  else 
  {  
    // Remove old map
    seq_data->model->shape_to_proxy->clear();

    // Add new proxies
    seq_data->model->structure->CreateNextLevelProxies( level_idx, 
                                               seq_data->parametrization.data(),
                                               optim_params->segment_length );

    // lets add features to proxies, so they cover approximate plane
    std::vector<Proxy*> top_level_proxies = optproxy::get_proxies_at_level(seq_data, 1);
    std::vector<R3Point> disk_points;
    for( r32 angle = 0 ; angle < 2 * RN_PI ; angle += RN_PI_OVER_TWO )
    {
      R3Point p( cos(angle), 0.0, sin(angle) );
      disk_points.push_back( p );
    }

    for ( Proxy * proxy : top_level_proxies )
    {
      r32 radius = 0.5;
      FETFeature * feature = proxy->feature;
      FETShape * shape     = feature->Shape();
      R3Point pos          = feature->Position(TRUE);
      R3Vector normal      = feature->Normal(TRUE);

      // Compute transformation 
      bsc::vec3d p( pos.X(), pos.Y(), pos.Z() );
      bsc::vec3d n( normal.X(), normal.Y(), normal.Z() );
      bsc::vec3d v = bsc::normalize( bsc::cross( n, bsc::vec3d( 0.0f, 1.0f, 0.0f ) ) );
      bsc::vec3d u = bsc::normalize( bsc::cross( v, n ) );
      bsc::mat4d xform;
      xform[0] = bsc::vec4d( u, 0.0f );
      xform[1] = bsc::vec4d( n, 0.0f );
      xform[2] = bsc::vec4d( v, 0.0f );
      xform[3] = bsc::vec4d( p, 1.0f );

      // Apply it to disk points and insert into proxy shapes
      R3Affine T ( bscmat4_to_R4Matrix(xform) );
      for ( i32 point_idx = 0; point_idx < disk_points.size() ; ++point_idx )
      {
        R3Point p = disk_points[point_idx] * radius;
        p.Transform( T );
        
        FETFeature * nf = new FETFeature ( *feature );
        nf->SetPosition( p );
        shape->InsertFeature ( nf );
      }
    }

    // update indices 
    for ( i32 i = 0; i < seq_data->model->structure->NProxies(); ++i )
    {
      Proxy * proxy = seq_data->model->structure->GetProxy( i );
      proxy->structure_index = i;
    }

    update_shape_proxy_map( seq_data->model );
  }
}


void
add_geometrical_correspondeces( const std::vector<StructRel> * relationships, 
                             std::vector<FETCorrespondence*> & corrs,
                             const Options * opts  )
{
  for ( int idx = 0;
        idx < relationships->size(); 
        idx++ )
  {
    StructRel sr = relationships->at( idx );

    FETShape *s1 = sr.proxy_A->shape;
    FETShape *s2 = sr.proxy_B->shape;
    i32 limit = bsc::min( s1->NFeatures(), s2->NFeatures() );

    bool inverse = false;
    if(sr.relationship_type == ANTIPARALLEL_RELATIONSHIP ||
       sr.relationship_type == PERPENDICULAR_RELATIONSHIP ) inverse = true;
    for( i32 i = 0 ; i < limit ; ++i )
    {
      FETFeature *f1 = NULL, *f2 = NULL;
      f1 = s1->Feature( i );
      int j = i;
      if( i > 0 && inverse ) j = limit - i;
      f2 = s2->Feature( j );
      if( f1 && f2 )
      {
        FETCorrespondence * corr = new FETCorrespondence( NULL,
                                                          f1,
                                                          f2,
                                                          sr.score,
                                                          sr.relationship_type);
        corrs.push_back( corr );
      }
    }
  }
  if( opts->print_verbose ) printf("      -- Added %lu geometrical"
                                   " correspondences\n", corrs.size() );
}


////////////////////////////////////////////////////////////////////////////////
// To Model Correspondences
////////////////////////////////////////////////////////////////////////////////

void random_shuffle( std::vector<int> & indices )
{
  size_t n = indices.size();
  for ( size_t i = n-1; i > 0; --i) 
  {
      std::swap( indices[i], indices[ get_random_int(0, 0, n-1) ] );
  }
}

void
add_hierarchical_correspondences( const FETReconstruction * reconstruction,
                                  const StructuralModel * model,
                                  std::vector<FETCorrespondence*> & corrs,
                                  const Options * opts )
{

  corrs.clear();
  for ( i32 ci = 0 ; ci < reconstruction->NCorrespondences() ; ci++ )
  {
    FETCorrespondence * corr = reconstruction->Correspondence( ci );

    FETFeature * f1 = corr->Feature(0);
    FETFeature * f2 = corr->Feature(1);

    auto sp1_it = model->shape_to_proxy->find( f1->Shape() );
    auto sp2_it = model->shape_to_proxy->find( f2->Shape() );
    
    i32 level1 = -1;
    i32 level2 = -1;
    if( sp1_it != model->shape_to_proxy->end() ) // just a feature in orginal reconstruction
    {
      Proxy *proxy = sp1_it->second;
      if( proxy ) level1 = proxy->level;
    }
    if( sp2_it != model->shape_to_proxy->end() ) // just a feature in original reconstruction
    {
      Proxy *proxy = sp2_it->second;
      if( proxy ) level2 = proxy->level;
    }

    // correspondences between higher level proxies
    if( level1 == 0 && level2 == 1  ) 
    {
      corrs.push_back( corr );
      FETShape * s1 = f1->Shape();
      FETShape * s2 = f2->Shape();
      assert( s1->NFeatures() < s2->NFeatures() );

      std::vector<int> indices( s1->NFeatures() );
      std::iota( indices.begin(), indices.end(), 0 );
      random_shuffle( indices );

      for (int fi = 0 ; fi < s2->NFeatures() ; ++fi )
      {
        FETFeature * f1 = s1->Feature( fi );
        FETFeature * f2 = s2->Feature( fi );
        FETCorrespondence * corr =  new FETCorrespondence( NULL, f1, f2, 1.0, COINCIDENT_RELATIONSHIP );
        corrs.push_back( corr );
      }
    }
    // correspondences between base level and shapes
    if( level1 == -1 && level2 == 0 )
    {
      R3Point pos = f1->Position(TRUE);
      R3Plane plane( f2->Position(TRUE), f2->Normal(TRUE) );
      pos.Project( plane );
      FETFeature * nf = new FETFeature( *f2 );
      f2->Shape()->InsertFeature( nf );
      nf->SetPosition( pos );
      FETCorrespondence * corr = 
            new FETCorrespondence( NULL, nf, f1, 1.0, COINCIDENT_RELATIONSHIP );
      corrs.push_back( corr );
    }
    if( level1 == -1 && level2 == -1 )
    {
      printf("This should not print\n");
    }
    if( level1 == 0 && level2 == -1 )
    {
      printf("This should not print\n");
    }
    if( level1 == 1 && level2 == 0 )
    {
      printf("This should not print\n");
    }
  }
  if( opts->print_verbose ) printf("      -- Added %lu hierarchical"
                                   " correspondences\n", corrs.size() );
}

}

////////////////////////////////////////////////////////////////////////////////
// Spatio-Temporal Processing
////////////////////////////////////////////////////////////////////////////////

r32 get_segment_extend( const bsc::mat4d * xforms,
                      i32 idx_A, i32 idx_B )
{
  r64 length = 0.0f;
  
  i32 min_idx = std::min( idx_A, idx_B );
  i32 max_idx = std::max( idx_A, idx_B );

  for ( i32 i = min_idx + 1; i <= max_idx; ++i )
  {
    bsc::mat4d xform_A = xforms[ i - 1 ];
    bsc::mat4d xform_B = xforms[ i ];
    bsc::vec3d pos_A = bsc::vec3d( xform_A[3] );
    bsc::vec3d pos_B = bsc::vec3d( xform_B[3] );
    length += bsc::norm( pos_A - pos_B );
  }

  return length;
}

void
get_segments_extends( std::vector<Segment> * segments,
                    const r32 max_dist_travelled,
                    const SequenceData * seq_data,
                    const Options * opts )
{
  using namespace bsc;

  // useful variables
  r32 current_length = 0.0f;
  r32 overlap_factor = 0.5f;
  segments->clear();

  // fill in the start location
  Segment s;
  s.start = 0;
  segments->push_back( s );
  for ( i32 i = 1; 
        i < opts->frame_end;
        ++i )
  {
    mat4d xform_A = seq_data->current_xforms[ i - 1 ];
    mat4d xform_B = seq_data->current_xforms[ i ];
    vec3d pos_A = vec3d( xform_A[3] );
    vec3d pos_B = vec3d( xform_B[3] );

    r32 length = norm( pos_A - pos_B );
    current_length += length;
    if( current_length > overlap_factor * max_dist_travelled )
    {
      Segment s;
      s.start = i;
      segments->push_back(s);
      current_length = 0.0f;
    }
  }

  // from each start location, count till reach required segment length
  for ( u32 i = 0 ;
        i < segments->size();
        ++i )
  {
    Segment * s = &(segments->at(i));
    s->length = 0.0f;
    for ( u32 j = s->start + 1;
          j < opts->frame_end;
          j++ )
    {
      mat4d xform_A = seq_data->current_xforms[ j - 1 ];
      mat4d xform_B = seq_data->current_xforms[ j ];
      vec3d pos_A = vec3d( xform_A[3] );
      vec3d pos_B = vec3d( xform_B[3] );
      r32 length = norm( pos_A - pos_B );

      s->length += length;
      if( s->length >= max_dist_travelled || j == opts->frame_end - 1 )
      {
        s->end = j;
        break;
      }
    }
  }

// deal with the last segment
  u32 n_segments = segments->size();
  if( segments->at( n_segments-1 ).end > opts->frame_end )
  {
    segments->at( n_segments-1  ).end = opts->frame_end - 1;
  }

  if( n_segments >= 2 )
  {
    if( segments->at( n_segments-1 ).length < 0.5f * max_dist_travelled )
    {
      Segment s = segments->at( n_segments-1 );
      segments->pop_back();

      Segment * t = &(segments->at( n_segments-2 ));
      t->end = s.end;
      t->length += s.length;
    }
  }
}

namespace opteq
{
////////////////////////////////////////////////////////////////////////////////
// Equations
////////////////////////////////////////////////////////////////////////////////


void 
add_trajectory_equations( FETReconstruction *reconstruction,
                        const SequenceData *data,
                        RNSystemOfEquations *system,
                        const OptimizationParameters *optim_params,
                        const Options *opts)
{
  // Check total weight
  if(optim_params->traj_weight <= 0)
    return;
  if(opts->frame_end == 0)
    return;
  if(optim_params->trajectory_sigma <= 0)
    return;

  // Get convenient variables
  RNScalar min_affinity = 1E-3;
  RNScalar f = -1.0 / (2.0 * optim_params->trajectory_sigma *
                             optim_params->trajectory_sigma);

  int extend = 12;

  // Determine total affinity
  RNScalar total_affinity = 0;
  float sigmoid_mean = 5.0f;
  float sigmoid_sigma = 0.5f;
  for (int i = 0; i < opts->frame_end; i++ )
  {
    r32 multiplier = mmath::sigmoid( data->pairwise_xforms_inliers[ i ], sigmoid_mean, sigmoid_sigma );
    for (int j = 1; j <= extend; j++ )
    {
      if(i + j < opts->frame_end)
      {
        RNScalar d = data->parametrization[i + j] - data->parametrization[i];
        RNScalar affinity = multiplier * exp(f * d * d);

        if(affinity < min_affinity) continue;
        total_affinity += affinity;
      }
    }
  }


  // Determine factor to compute weight per correspondence
  if(total_affinity > 0)
  {
    // Determine weighting scale factor
    RNScalar w = optim_params->traj_weight / total_affinity;

    // Add equations
    for (int i = 0; i < opts->frame_end; i++ )
    {
      r32 multiplier = mmath::sigmoid( data->pairwise_xforms_inliers[ i ], sigmoid_mean, sigmoid_sigma );
      FETShape *shape1 = reconstruction->Shape(i);
      R3Affine T1 ( bscmat4_to_R4Matrix( data->pairwise_xforms[i] ) );
      for (int j = 1; j <= extend; j++ )
      {
        if(i + j < opts->frame_end)
        {
          RNScalar d = data->parametrization[i + j] - data->parametrization[i];
          RNScalar affinity = multiplier * exp(f * d * d);
          if(affinity < min_affinity) continue;
          FETShape *shape2 = reconstruction->Shape(i + j);
          R3Affine T2 ( bscmat4_to_R4Matrix( data->pairwise_xforms[i+j] ) );

          R3Affine transformation21 = R3identity_affine;
          transformation21.InverseTransform(T1);
          transformation21.Transform(T2);
          reconstruction->AddPairwiseTransformationEquations( system, 
                                                              shape1, 
                                                              shape2, 
                                                              transformation21, 
                                                              w * affinity );
        }
      }
    }
  }
}

void
add_correspondence_equations( FETReconstruction *reconstruction,
                        RNSystemOfEquations *system,
                        const std::vector<FETCorrespondence *> &correspondences,
                        RNScalar total_weight,
                        const RNScalar * input_proportions,
                        std::string s, 
                        const Options *opts )
{
  // Check total weight
  if(total_weight <= 0)
    return;
  if(correspondences.size() == 0)
    return;

  RNScalar proportions[NUM_FEATURE_TYPES];
  if( input_proportions )
  {
    for (int i = 0; i < NUM_FEATURE_TYPES; ++i)
    {
      proportions[i] = input_proportions[i];
    }
  }
  else
  {
    // Uniform proportions.
    for (int i = 0; i < NUM_FEATURE_TYPES; ++i)
    {
      proportions[i] = 1.0f;
    }
  }

  // Compute total affinites
  RNScalar total_affinity = 0.0;
  RNScalar total_affinities[NUM_FEATURE_TYPES] = {0.0};
  for (int i = 0; i < correspondences.size(); i++)
  {
    FETCorrespondence *correspondence = correspondences[i];
    if( !correspondence ) continue;
    FETFeature *feature = correspondence->Feature(1);
    total_affinities[feature->GeneratorType()] += correspondence->affinity;
    total_affinity += correspondence->affinity;
  }

  for (int i = 0; i < NUM_FEATURE_TYPES; i++)
  {
    if(total_affinities[i] <= 0)
    {
      total_affinities[i] = 1.0;
      proportions[i] = 0.0;
    }
  }

  // normalize proportions
  RNScalar proportions_sum = 0.0;
  for (int i = 0; i < NUM_FEATURE_TYPES; i++)
    proportions_sum += proportions[i];
  for (int i = 0; i < NUM_FEATURE_TYPES; i++)
    proportions[i] /= proportions_sum;

  RNScalar total_weights[NUM_FEATURE_TYPES] = {0.0};
  for (int i = 0; i < NUM_FEATURE_TYPES; i++)
    total_weights[i] = proportions[i] * total_weight;

  // Determine factor to compute weight per correspondence
  RNScalar ws[NUM_FEATURE_TYPES] = {0.0};
  for (int i = 0; i < NUM_FEATURE_TYPES; i++)
    ws[i] = total_weights[i] / total_affinities[i];

  // Add equations for all correspondences
  RNScalar total_weights_check[NUM_FEATURE_TYPES] = {0.0};
  for (int i = 0; i < correspondences.size(); i++)
  {
    FETCorrespondence *correspondence = correspondences[i];
    if( !correspondence ) continue;
    FETFeature *feature1 = correspondence->Feature(0);
    FETFeature *feature2 = correspondence->Feature(1);
    RNScalar w = ws[feature2->GeneratorType()];
    if(w == 0)
    {
      continue;
    }
    if(!feature1 || !feature2)
      continue;
    FETShape *shape1 = feature1->shape;
    FETShape *shape2 = feature2->shape;
    if(!shape1 || !shape2)
      continue;
    RNScalar affinity = correspondence->affinity;
    if(w * affinity == 0)
    {
      printf("Error %f %f %d %d\n", w, affinity,
             feature1->GeneratorType(), feature2->GeneratorType());
      continue;
    }
    total_weights_check[feature2->GeneratorType()] += w * affinity;

    // Add equations for correspondence
    if(correspondence->relationship_type == COINCIDENT_RELATIONSHIP)
    {
      if((feature1->shape_type == POINT_FEATURE_SHAPE) &&
          (feature2->shape_type == POINT_FEATURE_SHAPE))
      {
        reconstruction->AddPointPointCorrespondenceEquations(system, shape1, shape2, feature1->Position(), feature2->Position(), w * affinity);
      }
      else
      {
        if(feature1->shape_type == LINE_FEATURE_SHAPE)
        {
          reconstruction->AddPointLineCorrespondenceEquations(system, shape2, shape1, feature2->Position(), feature1->Position(), feature1->Direction(), w * affinity);
        }
        else if(feature1->shape_type == PLANE_FEATURE_SHAPE)
        {
          reconstruction->AddPointPlaneCorrespondenceEquations(system, shape2, shape1, feature2->Position(), feature1->Position(), feature1->Normal(), w * affinity);
        }
        if(feature2->shape_type == LINE_FEATURE_SHAPE)
        {
          reconstruction->AddPointLineCorrespondenceEquations(system, shape1, shape2, feature1->Position(), feature2->Position(), feature2->Direction(), w * affinity);
        }
        else if(feature2->shape_type == PLANE_FEATURE_SHAPE)
        {
          reconstruction->AddPointPlaneCorrespondenceEquations(system, shape1, shape2, feature1->Position(), feature2->Position(), feature2->Normal(), w * affinity);
        }
      }
    }
    else if(correspondence->relationship_type == PARALLEL_RELATIONSHIP)
    {
      if(feature1->shape_type == LINE_FEATURE_SHAPE)
      {
        if(feature2->shape_type == LINE_FEATURE_SHAPE)
        {
          reconstruction->AddParallelVectorEquations(system, shape1, shape2, feature1->Direction(), feature2->Direction(), w * affinity);
        }
        else if(feature2->shape_type == PLANE_FEATURE_SHAPE)
        {
          reconstruction->AddPerpendicularVectorEquations(system, shape1, shape2, feature1->Direction(), feature2->Normal(), w * affinity);
        }
      }
      else if(feature1->shape_type == PLANE_FEATURE_SHAPE)
      {
        if(feature2->shape_type == LINE_FEATURE_SHAPE)
        {
          reconstruction->AddPerpendicularVectorEquations(system, shape1, shape2, feature1->Normal(), feature2->Direction(), w * affinity);
        }
        else if(feature2->shape_type == PLANE_FEATURE_SHAPE)
        {
          reconstruction->AddParallelVectorEquations(system, shape1, shape2, feature1->Normal(), feature2->Normal(), w * affinity);
        }
      }
    }
    else if(correspondence->relationship_type == ANTIPARALLEL_RELATIONSHIP)
    {
      if(feature1->shape_type == LINE_FEATURE_SHAPE)
      {
        if(feature2->shape_type == LINE_FEATURE_SHAPE)
        {
          reconstruction->AddParallelVectorEquations(system, shape1, shape2, feature1->Direction(), -feature2->Direction(), w * affinity);
        }
        else if(feature2->shape_type == PLANE_FEATURE_SHAPE)
        {
          reconstruction->AddPerpendicularVectorEquations(system, shape1, shape2, feature1->Direction(), feature2->Normal(), w * affinity);
        }
      }
      else if(feature1->shape_type == PLANE_FEATURE_SHAPE)
      {
        if(feature2->shape_type == LINE_FEATURE_SHAPE)
        {
          reconstruction->AddPerpendicularVectorEquations(system, shape1, shape2, feature1->Normal(), feature2->Direction(), w * affinity);
        }
        else if(feature2->shape_type == PLANE_FEATURE_SHAPE)
        {
          reconstruction->AddParallelVectorEquations(system, shape1, shape2, feature1->Normal(), -feature2->Normal(), w * affinity);
        }
      }
    }
    else if(correspondence->relationship_type == PERPENDICULAR_RELATIONSHIP)
    {
      if(feature1->shape_type == LINE_FEATURE_SHAPE)
      {
        if(feature2->shape_type == LINE_FEATURE_SHAPE)
        {
          reconstruction->AddPerpendicularVectorEquations(system, shape1, shape2, feature1->Direction(), feature2->Direction(), w * affinity);
        }
        else if(feature2->shape_type == PLANE_FEATURE_SHAPE)
        {
          reconstruction->AddParallelVectorEquations(system, shape1, shape2, feature1->Direction(), feature2->Normal(), w * affinity);
        }
      }
      else if(feature1->shape_type == PLANE_FEATURE_SHAPE)
      {
        if(feature2->shape_type == LINE_FEATURE_SHAPE)
        {
          reconstruction->AddParallelVectorEquations(system, shape1, shape2, feature1->Normal(), feature2->Direction(), w * affinity);
        }
        else if(feature2->shape_type == PLANE_FEATURE_SHAPE)
        {
          reconstruction->AddPerpendicularVectorEquations(system, shape1, shape2, feature1->Normal(), feature2->Normal(), w * affinity);
        }
      }
    }
  }
  if( opts->print_verbose )
  {
    printf("     Total Weights for %15s Correspondences : "
          "%11.3f %11.3f %11.3f %11.3f %11.3f "
          "%11.3f %11.3f %11.3f %11.3f %11.3f",
          s.c_str(),
          total_weights_check[SIFT_FEATURE_TYPE],
          total_weights_check[FAST_FEATURE_TYPE],
          total_weights_check[CORNER_FEATURE_TYPE],
          total_weights_check[RIDGE_FEATURE_TYPE],
          total_weights_check[VALLEY_FEATURE_TYPE],
          total_weights_check[SILHOUETTE_FEATURE_TYPE],
          total_weights_check[BORDER_FEATURE_TYPE],
          total_weights_check[UNIFORM_FEATURE_TYPE],
          total_weights_check[PLANE_FEATURE_TYPE],
          total_weights_check[STRUCTURE_FEATURE_TYPE]);
  }
}
}

////////////////////////////////////////////////////////////////////////////////
// Main Optimization
////////////////////////////////////////////////////////////////////////////////

void
optimize_transformations(
    FETReconstruction *reconstruction,
    const SequenceData *seq_data,
    const OptimizationParameters *optim_params,
    const Options *opts)
{

  // Check the shapes
  if(reconstruction->NShapes() < 2)
    return;

  // Update variable indices
  int n = 0;
  RNTime t;

  for (int i = 0; i < reconstruction->NShapes(); i++)
  {
    FETShape *shape = reconstruction->Shape(i);
    shape->UpdateVariableIndex(n);
  }

  // Create system of equations
  RNSystemOfEquations *equations = new RNSystemOfEquations(n);
  
  reconstruction->AddInertiaEquations(equations,
                                      optim_params->inertia_weight);

  if(optim_params->traj_weight > 0)
  {
    opteq::add_trajectory_equations(reconstruction,
                           seq_data,
                           equations,
                           optim_params,
                           opts);
  }

  if( !seq_data->cp_corrs.empty() || 
       !seq_data->hierarchical_corrs.empty() || 
       !seq_data->geometrical_corrs.empty() )
  {    
    if( opts->print_verbose )
    {
      printf("                                                         "
            "%11s %11s %11s %11s %11s %11s %11s %11s %11s %11s\n",
            "Sift",
            "Fast",
            "Corner",
            "Ridge",
            "Valley",
            "Silhouette",
            "Border",
            "Uniform",
            "Plane",
            "Structure" );
    }
  }

  if( optim_params->cp_weight > 0.0f && 
       !seq_data->cp_corrs.empty() )
  {
    t.Read();
    opteq::add_correspondence_equations(reconstruction,
                                        equations,
                                        seq_data->cp_corrs,
                                        optim_params->cp_weight,
                                        optim_params->cp_proportions,
                                        " Closest Point ",
                                        opts );
    if( opts->print_verbose ) printf("| Time : %f\n", t.Elapsed() );
  }

  if( optim_params->hierarchical_weight > 0.0f && 
       !seq_data->hierarchical_corrs.empty() )
  {        
    t.Read();
    opteq::add_correspondence_equations(reconstruction,
                                        equations,
                                        seq_data->hierarchical_corrs,
                                        optim_params->hierarchical_weight, 
                                        NULL,
                                        " Hierarchy ", 
                                        opts);
    if( opts->print_verbose ) printf("| Time : %f\n", t.Elapsed() );        }

if( optim_params->geometrical_weight > 0.0f && 
       !seq_data->geometrical_corrs.empty() )
  {
    t.Read();
    opteq::add_correspondence_equations(reconstruction,
                                        equations,
                                        seq_data->geometrical_corrs,
                                        optim_params->geometrical_weight, 
                                        NULL,   
                                        " Geometry ", 
                                        opts );
    if( opts->print_verbose ) printf("| Time : %f\n", t.Elapsed() );
  }

  // Initialize variables
  double *x = new double[n];
  for (int i = 0; i < n; i++)
    x[i] = 0;
  if( opts->print_verbose )
  {
    printf("     NVars: %d | NEqn: %d | NParDeriv: %d -> Residuals before: %g ",
            equations->NVariables(), 
            equations->NEquations(), 
            equations->NPartialDerivatives(), 
            sqrt(equations->SumOfSquaredResiduals(x)));
  }

  // Solve system of equations
  if(equations->NEquations() >= n)
  {
    RNTime t_minimization; t_minimization.Read();
    if(!equations->Minimize(x, reconstruction->solver))
    {
      fprintf(stderr, "Unable to minimize system of equations\n");
      delete[] x;
      return;
    }
    if( opts->print_verbose )
    {
      printf("after: %g | Time to minimize system of equations: %f \n", 
          sqrt(equations->SumOfSquaredResiduals(x)), t_minimization.Elapsed() );
    }
  }

  // Extract solution
  for (int i = 0; i < reconstruction->NShapes(); i++)
  {
    FETShape *shape = reconstruction->Shape(i);
    shape->UpdateVariableValues(x);
  }

  // Delete variables
  delete[] x;
  delete equations; // <- This might be calling a TON of destructors
}

#endif
