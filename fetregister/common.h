// FineToCoarseRegistration - fetregister - common file. 
//
// This file includes definitions of main datatypes used by the system, as well
// as some random helper functions for:
// - argument parsing
// - calculating current parametrization length
// - realigning registration with main world axes
// - distributions for random number generations


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
// Inludes
////////////////////////////////////////////////////////////////////////////////


#define BSC_IMPLEMENTATION
#ifdef FETREGISTER_USE_WINDOW
  #define BSC_USE_WINDOW
  #define BSC_USE_IMGUI
#endif
#include "basics.h"

#include "FET/FET.h"

#include <random>
#include <numeric>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <thread>
#include <set>
#include <cstdlib>

#include "structure.h"

////////////////////////////////////////////////////////////////////////////////
// Global consts
////////////////////////////////////////////////////////////////////////////////

// TODO(maciej): Possibly get rid of those and replace with std::vector?
#define MAX_STR_LENGTH        512
#define MAX_FRAME_NAME_LENGTH 64
#define MAX_ITER              32
#define MAX_N_GT_CORRS        2048

////////////////////////////////////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////////////////////////////////////

enum FeatureDisplayMode
{
  SHOW_COLORS = 0,
  SHOW_GENERATOR_TYPE = 1,
  SHOW_PRIMITIVE_IDX = 2,
  SHOW_NORMALS = 3,
  SHOW_PHONG = 4,
};

enum ProxyDisplayMode
{
  SHOW_LEVEL = 0,
  SHOW_PARENT_IDX = 1,
  SHOW_PROXY_IDX = 2
};

struct Correspondence
{
  char *name_A, *name_B;
  i32 idx_A, idx_B;
  bsc::vec3d point_A, point_B;
};

struct Segment
{
  i32 start;
  i32 end;
  r32 length;
};

struct SeenInfo
{
  int count;
  bool incremented;
};

struct StructRel
{
  Proxy * proxy_A;
  Proxy * proxy_B;
  int relationship_type;
  float score;
};

struct StructuralModel
{
  Structure *structure = NULL;
  std::unordered_map<FETShape*, Proxy*> *shape_to_proxy;
  std::vector<StructRel> *relationships;
};

struct Options
{
  char * conf_filename;
  char * fet_filename;
  char * output_conf_filename;
  char * output_error_filename;
  char * registration_conf_filename;

  i32 frame_end;
  i32 downsample_factor;
  i32 solver;

  bool print_verbose;
  bool load_frames;
  bool batch_mode;
  bool screenshot_mode;
  bool write_ply_model;
};

struct ViewOptions
{
  bool show_frames                  = 1;
  bool show_cp_corrs                = 0;
  bool show_gt_corrs                = 0;
  bool show_geometrical_corrs       = 0;
  bool show_hierarchical_corrs      = 0;
  bool show_cameras                 = 0;
  bool show_bboxes                  = 0;
  bool show_proxies                 = 1;
  bool show_matches                 = 0;
  bool show_trajectory              = 0;
  bool show_grid                    = 1;
  bool show_axes                    = 1;
  bool show_fps                     = 1;
  bool inspection_mode              = 0;
  int selected_class_id             = -1;
  int selected_object_id            = -1;

  r32 fov_y                         = 20.0f;
  r32 point_size                    = 6.0f;
  r32 line_width                    = 2.0f;
  r32 max_dist                      = 5.0f;
  r32 max_height                    = 3.0f;

  i32 frame_idx_A                   = 0;
  i32 frame_idx_B                   = 0;
  i32 frame_skip                    = 0;
  bool show_only_pair               = 0;

  bsc::mat4d * xforms_ptr           = NULL;
  i32 xform_idx                     = 0;

  FeatureDisplayMode feat_disp_mode = SHOW_COLORS;
  ProxyDisplayMode proxy_disp_mode  = SHOW_LEVEL;
  
  r32 proxy_offset                  = 0.0f;
  i32 proxy_level                   = -1;
  i32 xform_scheme                  = CURRENT_TRANSFORMATION;

  bool show_relationships_by_type[4];
};

struct OptimizationParameters
{
  r32 trajectory_sigma;
  r32 angle_threshold;
  r32 distance_threshold;
  r32 init_angle_threshold;
  r32 final_angle_threshold;
  r32 init_distance_threshold;
  r32 final_distance_threshold;

  r32 structure_angle_threshold;
  r32 structure_distance_threshold;
  
  r32 init_cp_weight;
  r32 init_traj_weight;
  r32 init_inertia_weight;
  r32 init_geometrical_weight;
  r32 init_hierarchical_weight;

  r32 final_cp_weight;
  r32 final_traj_weight;
  r32 final_inertia_weight;
  r32 final_geometrical_weight;
  r32 final_hierarchical_weight;

  r32 cp_weight;
  r32 traj_weight;
  r32 inertia_weight;
  r32 geometrical_weight;
  r32 hierarchical_weight;

  i32 n_samples;
  i32 max_n_corrs;
  r32 max_depth;

  i32 n_iter;
  i32 iter_idx;
  r32 segment_length;
  r32 segment_length_growth;

  i32 n_threads;

  r32 glob_max_dist, glob_min_dist;
  r32 glob_max_angle, glob_min_angle;

  r64 cp_proportions[ NUM_FEATURE_TYPES ];
  std::unordered_map<i32, SeenInfo> seen_count;
};

static std::vector<std::mt19937> gens;
static std::vector<int> counters;

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
  bsc::mat4d * rgbdsfm_xforms               = NULL;
  bsc::mat4d * robust_reconstruction_xforms = NULL;
  bsc::mat4d * elastic_fusion_xforms        = NULL;
  bsc::mat4d * kintinuous_xforms            = NULL;
  bsc::mat4d * gt_xforms                    = NULL;

  bsc::mat4d * xforms_history               = NULL;
  i32 n_xforms                              = 0;

  // reconstruction stuff
  FETReconstruction * reconstruction = NULL;
  StructuralModel * model            = NULL;
  R3Box * bboxes                     = NULL;

  // correspondence storage
  std::vector<FETCorrespondence*> cp_corrs;
  std::vector<FETCorrespondence*> geometrical_corrs;
  std::vector<FETCorrespondence*> hierarchical_corrs;

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


void 
initialize_options( Options * opts,
                    ViewOptions * view_opts,
                    OptimizationParameters * optim_params )
{
  opts->conf_filename                                     = NULL;
  opts->output_conf_filename                              = NULL;
  opts->output_error_filename                             = NULL;
  opts->registration_conf_filename                        = NULL;

  opts->frame_end                                         = 0;

  opts->downsample_factor                                 = 2;
  opts->solver                                            = RN_CSPARSE_SOLVER;

  opts->print_verbose                                     = false;
  opts->load_frames                                       = false;
  opts->batch_mode                                        = false;
  opts->screenshot_mode                                   = false;
  opts->write_ply_model                                   = false;

  view_opts->show_frames                                  = 1;
  view_opts->show_cp_corrs                                = 0;
  view_opts->show_gt_corrs                                = 0;
  view_opts->show_geometrical_corrs                       = 0;
  view_opts->show_hierarchical_corrs                      = 0;
  view_opts->show_cameras                                 = 0;
  view_opts->show_bboxes                                  = 0;
  view_opts->show_proxies                                 = 0;
  view_opts->show_matches                                 = 0;
  view_opts->show_trajectory                              = 1;
  view_opts->show_grid                                    = 1;
  view_opts->show_axes                                    = 0;
  view_opts->show_fps                                     = 1;
  view_opts->inspection_mode                              = 0;
  view_opts->selected_class_id                            = -1;
  view_opts->selected_object_id                           = -1;

  view_opts->fov_y                                        = 45.0f;
  view_opts->point_size                                   = 6.0f;
  view_opts->line_width                                   = 2.0f;
  view_opts->max_dist                                     = 3.5f;
  view_opts->max_height                                   = 3.0f;

  view_opts->frame_idx_A                                  = 0;
  view_opts->frame_idx_B                                  = 0;
  view_opts->show_only_pair                               = 0;

  view_opts->xforms_ptr                                   = NULL;
  view_opts->xform_idx                                    = 0;

  view_opts->feat_disp_mode                               = SHOW_COLORS;
  view_opts->proxy_disp_mode                              = SHOW_LEVEL;
  view_opts->proxy_level                                  = -1;
  view_opts->proxy_offset                                 = 0.0f;
  view_opts->xform_scheme                                 = CURRENT_TRANSFORMATION;

  view_opts->show_relationships_by_type[0]                = 1;
  view_opts->show_relationships_by_type[1]                = 1;
  view_opts->show_relationships_by_type[2]                = 1;
  view_opts->show_relationships_by_type[3]                = 1;

  optim_params->trajectory_sigma                          = 0.25;
  optim_params->angle_threshold                           = 0.0;
  optim_params->distance_threshold                        = 0.2;
  optim_params->init_angle_threshold                      = 0.5;
  optim_params->final_angle_threshold                     = 0.35;
  optim_params->init_distance_threshold                   = 0.5;
  optim_params->final_distance_threshold                  = 0.2;

  optim_params->structure_angle_threshold                 = 0.15 * RN_PI;
  optim_params->structure_distance_threshold              = 0.15;

  optim_params->init_cp_weight                            = 1500.0f;
  optim_params->init_traj_weight                          = 1000.0f;
  optim_params->init_geometrical_weight                   = 1500.0f;
  optim_params->init_hierarchical_weight                  = 1500.0f;
  optim_params->init_inertia_weight                       = 1.0f;

  optim_params->final_cp_weight                           = 1000.0f;
  optim_params->final_traj_weight                         = 1000.0f;
  optim_params->final_geometrical_weight                  = 1000.0f;
  optim_params->final_hierarchical_weight                 = 1000.0f;
  optim_params->final_inertia_weight                      = 1.0f;

  optim_params->cp_weight                                 = -1.0f;
  optim_params->traj_weight                               = -1.0f;
  optim_params->geometrical_weight                        = -1.0f;
  optim_params->inertia_weight                            = -1.0f;

  optim_params->max_n_corrs                               = 64000;
  optim_params->n_samples                                 = 256;
  optim_params->max_depth                                 = 5.0f;

  optim_params->n_iter                                    = 12;
  optim_params->iter_idx                                  = 0;
  optim_params->segment_length                            = 3.0f;
  optim_params->segment_length_growth                     = 2.0f;

  optim_params->cp_proportions[ UNKNOWN_FEATURE_TYPE ]    = 0.0;
  optim_params->cp_proportions[ BORDER_FEATURE_TYPE ]     = 0.0;
  optim_params->cp_proportions[ SHADOW_FEATURE_TYPE ]     = 0.0;
  optim_params->cp_proportions[ POLE_FEATURE_TYPE ]       = 0.0;
  optim_params->cp_proportions[ SIFT_FEATURE_TYPE ]       = 0.0;
  optim_params->cp_proportions[ FAST_FEATURE_TYPE ]       = 0.0;
  optim_params->cp_proportions[ CORNER_FEATURE_TYPE ]     = 0.0;
  optim_params->cp_proportions[ RIDGE_FEATURE_TYPE ]      = 0.2;
  optim_params->cp_proportions[ VALLEY_FEATURE_TYPE ]     = 0.2;
  optim_params->cp_proportions[ SILHOUETTE_FEATURE_TYPE ] = 0.6;
  optim_params->cp_proportions[ UNIFORM_FEATURE_TYPE ]    = 0.0;
  optim_params->cp_proportions[ PLANE_FEATURE_TYPE ]      = 1.0;
  optim_params->cp_proportions[ STRUCTURE_FEATURE_TYPE  ] = 0.0;

  optim_params->n_threads = 1;
};

void parse_arguments( int argc, char **argv, 
                     Options *opts,
                     OptimizationParameters *optim_params )
{
  bsc::arg_parse args;

  args.add(bsc::argument<char *>("fet_filename",
                                 "Input reconstruction filename",
                                 &(opts->fet_filename)));

  args.add(bsc::argument<char *>("conf_filename",
                                 "Input configuration filename",
                                 &(opts->conf_filename)));

  args.add(bsc::argument<bool>("-v",
                               "--verbose",
                               "Prints additional information",
                               &(opts->print_verbose), 0));

  args.add(bsc::argument<bool>("-p",
                               "--output_ply",
                               "Writes ply file to disk. Name based on "
                               "optimization params",
                               &(opts->write_ply_model), 0));

  args.add(bsc::argument<bool>("-b",
                               "--batch",
                               "Run in a batch mode",
                               &(opts->batch_mode), 0));

  args.add(bsc::argument<bool>("-t",
                               "--screenshot",
                               "Run in a screenshot mode",
                               &(opts->screenshot_mode), 0));


  args.add(bsc::argument<r32>( "-iw",
                               "--initial_weights",
                               "Initial weights (Closest Points, Trajectory, "
                               " Inertia, Hierarchy, Geometry)",
                                &(optim_params->init_cp_weight), 5 ) );

  args.add(bsc::argument<r32>( "-fw",
                               "--final_weights",
                               "Final weights (Closest Points, Trajectory, "
                               " Inertia, Hierarchy, Geometry)",
                                &(optim_params->final_cp_weight), 5 ) );

  args.add(bsc::argument<r32>( "-md",
                               "--max_depth",
                               "Maximum allowed distance from camera",
                               &(optim_params->max_depth) ) );

  args.add(bsc::argument<r32>( "-at",
                               "--angle_threshold",
                               "Initial and final angle threshold in redians",
                               &(optim_params->init_angle_threshold), 2));

  args.add(bsc::argument<r32>( "-dt",
                               "--distance_threshold",
                               "Initial and final distance threshold in meters",
                               &(optim_params->init_distance_threshold), 2));

  args.add(bsc::argument<r32>( "-sl",
                               "--segment_length",
                               "Initial segment length in meters",
                               &(optim_params->segment_length) ) );

  args.add(bsc::argument<r32>( "-sg",
                               "--segment_length_growth",
                               "How much should segment length increase from "
                               "iteration to iteration",
                               &(optim_params->segment_length_growth) ) );

  args.add(bsc::argument<i32>( "-ni",
                               "--num_iterations",
                               "Number of iterations the optimization will "
                               "be run",
                               &(optim_params->n_iter) ) );

  args.add(bsc::argument<i32>( "-cs",
                               "--n_corr_samples",
                               "Maximum number of correspondence samples",
                               &(optim_params->n_samples) ) );

  args.add(bsc::argument<i32>("-ef",
                              "--end_frame",
                              "Index of the last frame to load",
                              &(opts->frame_end)));

  args.add(bsc::argument<char *>("-o",
                                 "--output_conf",
                                 "Filename of output .conf file",
                                 &(opts->output_conf_filename)));

  args.add(bsc::argument<int>( "-nt",
                               "--n_threads",
                               "Number of threads",
                               &(optim_params->n_threads)));

  args.add(bsc::argument<char *>("-e",
                                 "--output_errors",
                                 "Filename of error file",
                                 &(opts->output_error_filename)));

  args.parse(argc, argv);
}

void
initialize_random_numbers( OptimizationParameters * optim_params )
{
  unsigned long long base_seed = 12346ULL;
  for ( int i = 0; i < optim_params->n_threads; ++i )
  {
    gens.push_back( std::mt19937( base_seed + i ) );
  }
}


namespace rand_dist
{
// This is from libc++ and http://www.blahonga.org/2015/04/state-of-art-of-random-double-generation.html
// Goal is to ensure that the same sequences of random doubles are generated on
// all the systems. This is a replacement for std::uniform_real_distribution,
// which produce different sequences on UNIX/MACOSX.
// Moreover given that we only want numbers in range 0 to 1 we can simplify this
// However, the link above suggest that std::uniform_real_distribution is flawed.
// TODO(maciej): Investigate the alternate uniform rands.

double
generate_canonical(std::mt19937 &gen)
{
  // const size_t Dt = std::numeric_limits<double>::digits;
  // const size_t b = Dt < bits ? Dt : bits;
  // const size_t logR = std::log2( gen.max() + uint64_t(1) );
  const size_t k = 2; // b / logR + (b % logR != 0) + (b == 0);
  const double Rp = gen.max() + double(1);
  double base = Rp;
  double Sp = gen();
  for (size_t i = 1; i < k; ++i, base *= Rp)
      Sp += gen() * base;
  return Sp / base;
}

// Below code is also from clangs libc++ - it is simplified version for only
// std::mt19937 as random engine and int as a return type(removed templating)
// It is bit crazy but we want uniform_int_distribution with a sequence that is
// replicated across compilers.

#ifndef _MSC_VER
uint32_t inline
clz(uint32_t x)
{
    return static_cast<uint32_t>(__builtin_clz(x));
}
#else
//NOTE(maciej): I have not tested msvc, just gcc and clang.
//This is from http://stackoverflow.com/questions/355967/how-to-use-msvc-intrinsics-to-get-the-equivalent-of-this-gcc-code
#include <intrin.h>
uint32_t __inline 
clz( uint32_t value )
{
    DWORD leading_zero = 0;

    if ( _BitScanReverse( &leading_zero, value ) )
    {
      return static_cast<uint32_t>(31 - leading_zero);
    }
    else
    {
      // Same remarks as above
      return static_cast<uint32_t>(32);
    }
}
#endif

// independent_bits_engine
class independent_bits_engine
{
private:
    std::mt19937& e_;
    size_t w_;
    size_t w0_;
    size_t n_;
    size_t n0_;
    uint32_t y0_;
    uint32_t y1_;
    uint32_t mask0_;
    uint32_t mask1_;

    static const uint32_t Rp = 0; //std::mt19937::max() - std::mt19937::min() + uint32_t(1);
    static const size_t m    = 32;//log2<uint32_t, _Rp>::value;
    static const size_t WDt  = 32;//std::numeric_limits<uint32_t>::digits;
    static const size_t EDt  = 32;//std::numeric_limits<uint32_t>::digits;

public:
    // constructors and seeding functions
    independent_bits_engine(std::mt19937& e, size_t w);

    // generating functions
    uint32_t operator()() {return eval(std::integral_constant<bool, Rp != 0>());}

private:
    uint32_t eval(std::false_type);
    uint32_t eval(std::true_type);
};

independent_bits_engine
    ::independent_bits_engine(std::mt19937& e, size_t w)
        : e_(e),
          w_(w)
{
  n_ = w_ / m + (w_ % m != 0);
  w0_ = w_ / n_;
  if (Rp == 0) { y0_ = Rp; }
  else if (w0_ < WDt) { y0_ = (Rp >> w0_) << w0_; }
  else { y0_ = 0; }

  if (Rp - y0_ > y0_ / n_)
  {
      ++n_;
      w0_ = w_ / n_;
      if (w0_ < WDt)
          y0_ = (Rp >> w0_) << w0_;
      else
          y0_ = 0;
  }
  n0_ = n_ - w_ % n_;
  if (w0_ < WDt - 1) { y1_ = (Rp >> (w0_ + 1)) << (w0_ + 1); }
  else { y1_ = 0; }
  mask0_ = w0_ > 0 ? uint32_t(~0) >> (EDt - w0_) :
                      uint32_t(0);
  mask1_ = w0_ < EDt - 1 ? uint32_t(~0) >> (EDt - (w0_ + 1)) :
                            uint32_t(~0);
}

inline
uint32_t
independent_bits_engine::eval(std::false_type)
{
  return static_cast<uint32_t>(e_() & mask0_);
}

uint32_t
independent_bits_engine::eval(std::true_type)
{
  uint32_t Sp = 0;
  for (size_t k = 0; k < n0_; ++k)
  {
    uint32_t u;
    do
    {
      u = e_() - std::mt19937::min();
    } while (u >= y0_);
    if (w0_ < WDt) { Sp <<= w0_; }
    else { Sp = 0; }
    Sp += u & mask0_;
  }
  for (size_t k = n0_; k < n_; ++k)
  {
    uint32_t u;
    do
    {
      u = e_() - std::mt19937::min();
    } while (u >= y1_);
    if (w0_ < WDt - 1) { Sp <<= w0_ + 1; }
    else { Sp = 0; }
    Sp += u & mask1_;
  }
  return Sp;
}

// uniform_int_distribution

class uniform_int_distribution
{
public:
    // types

  class param_type
  {
    int a_;
    int b_;

    public:
      typedef uniform_int_distribution distribution_type;

      explicit param_type(int a = 0, int b = std::numeric_limits<int>::max())
          : a_(a), b_(b) {}

      int a() const {return a_;}
      int b() const {return b_;}

      friend bool operator==(const param_type& x, const param_type& y)
          {return x.a_ == y.a_ && x.b_ == y.b_;}
      friend bool operator!=(const param_type& x, const param_type& y)
          {return !(x == y);}
  };

private:
    param_type p_;

public:
    // constructors and reset functions
    explicit uniform_int_distribution(int a = 0, int b = std::numeric_limits<int>::max())
        : p_(param_type(a, b)) {}
    explicit uniform_int_distribution(const param_type& p) : p_(p) {}
    void reset() {}

    // generating functions
    int operator()(std::mt19937& g)
        {return (*this)(g, p_);}
    int operator()(std::mt19937& g, const param_type& p);

    // property functions
    int a() const {return p_.a();}
    int b() const {return p_.b();}

    param_type param() const {return p_;}
    void param(const param_type& p) {p_ = p;}

    int min() const {return a();}
    int max() const {return b();}

    friend bool operator==(const uniform_int_distribution& x,
                           const uniform_int_distribution& y)
        {return x.p_ == y.p_;}
    friend bool operator!=(const uniform_int_distribution& x,
                           const uniform_int_distribution& y)
            {return !(x == y);}
};

int
uniform_int_distribution::operator()(std::mt19937& g, const param_type& p)
{
    const uint32_t Rp = p.b() - p.a() + uint32_t(1);
    if (Rp == 1) { return p.a(); }
    const size_t Dt = 32;//std::numeric_limits<uint32_t>::digits;
    if (Rp == 0) { return static_cast<int>(independent_bits_engine(g, Dt)()); }
    size_t w = Dt - clz(Rp) - 1;
    if ((Rp & (std::numeric_limits<uint32_t>::max() >> (Dt - w))) != 0)
      ++w;
    independent_bits_engine e(g, w);
    uint32_t u;
    do
    {
        u = e();
    } while (u >= Rp);
    return static_cast<int>(u + p.a());
}

} //end rand_dist



double get_random_real( int thread_idx, double min = 0.0, double max = 1.0 )
{
  double retval = (max - min) * rand_dist::generate_canonical( gens[thread_idx] ) + min;
  return retval;
}

int get_random_int( int thread_idx, int min, int max )
{
  rand_dist::uniform_int_distribution idist(min, max);
  return idist(gens[thread_idx]);
}


////////////////////////////////////////////////////////////////////////////////
// Type Conversions
////////////////////////////////////////////////////////////////////////////////

static inline R4Matrix 
bscmat4_to_R4Matrix( const bsc::mat4d & im )
{
  R4Matrix om;
  om[0][0]=im[0].x; om[1][0]=im[0].y; om[2][0]=im[0].z; om[3][0]=im[0].w;
  om[0][1]=im[1].x; om[1][1]=im[1].y; om[2][1]=im[1].z; om[3][1]=im[1].w;
  om[0][2]=im[2].x; om[1][2]=im[2].y; om[2][2]=im[2].z; om[3][2]=im[2].w;
  om[0][3]=im[3].x; om[1][3]=im[3].y; om[2][3]=im[3].z; om[3][3]=im[3].w;
  return om;
}

static inline bsc::mat4d
R4Matrix_to_bscmat4( const R4Matrix & im )
{
  bsc::mat4d om;
  om[0].x=im[0][0]; om[0].y=im[1][0]; om[0].z=im[2][0]; om[0].w=im[3][0];
  om[1].x=im[0][1]; om[1].y=im[1][1]; om[1].z=im[2][1]; om[1].w=im[3][1];
  om[2].x=im[0][2]; om[2].y=im[1][2]; om[2].z=im[2][2]; om[2].w=im[3][2];
  om[3].x=im[0][3]; om[3].y=im[1][3]; om[3].z=im[2][3]; om[3].w=im[3][3];
  return om;
}

static inline bsc::vec3d
R3Point_to_bscvec3( const R3Point & p )
{
  return bsc::vec3d( p.X(), p.Y(), p.Z() );
}

////////////////////////////////////////////////////////////////////////////////
// Comparators for qsort
////////////////////////////////////////////////////////////////////////////////

int std_compare (const void * a, const void * b)
{
  return  (*(int*)a > *(int*)b) - (*(int*)a < *(int*)b);
}

int std_comparef (const void * a, const void * b)
{
  return  (*(r32*)a > *(r32*)b) - (*(r32*)a < *(r32*)b);
}

int rev_compare (const void * a, const void * b)
{
  return  (*(int*)a < *(int*)b) - (*(int*)a > *(int*)b);
}

struct tuple { float val; int idx; };
int idx_compare (const void * a, const void * b)
{
  tuple * t_a = (tuple*)a;
  tuple * t_b = (tuple*)b;

  return  (t_a->val < t_b->val) - (t_a->val > t_b->val);
}

// Taken from stack overflow
// http://stackoverflow.com/questions/24586499/keeping-track-of-the-original-ind-of-an-array-after-sorting-in-c
//-------------------------------
// maybe it is better to do it with c++11 features.
r32 *array_helper = NULL; // it needs this in outher scope
int arr_idx_compare( const void *a, const void *b ){
    int ia = *(int *)a;
    int ib = *(int *)b;
    return (array_helper[ia] < array_helper[ib]) ? 
                                    -1 : (array_helper[ia] > array_helper[ib]);
}
//-------------------------

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
      printf("common.h | calculate_rmse : This should never print!\n");
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
                       r32 *rmse, r32 *mean, r32 *sdev )
{
  rmse[0] = calculate_rmse( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                            seq_data->current_xforms );
  mean[0] = calculate_mean( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                            seq_data->current_xforms );
  sdev[0] = calculate_stddev( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              mean[0], seq_data->current_xforms );

  rmse[1] = NAN; mean[1] = NAN; sdev[1] = NAN;
  if ( seq_data->rgbdsfm_xforms )
  {
    rmse[1] = calculate_rmse( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              seq_data->rgbdsfm_xforms );
    mean[1] = calculate_mean( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              seq_data->rgbdsfm_xforms );
    sdev[1] = calculate_stddev( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                                mean[1],seq_data->rgbdsfm_xforms );
  }

  rmse[2] = NAN; mean[2] = NAN; sdev[2] = NAN;
  if ( seq_data->robust_reconstruction_xforms )
  {
    rmse[2] = calculate_rmse( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              seq_data->robust_reconstruction_xforms );
    mean[2] = calculate_mean( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              seq_data->robust_reconstruction_xforms );
    sdev[2] = calculate_stddev( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                            mean[2], seq_data->robust_reconstruction_xforms );
  }

  rmse[3] = NAN; mean[3] = NAN; sdev[3] = NAN;
  if ( seq_data->elastic_fusion_xforms )
  {
    rmse[3] = calculate_rmse( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              seq_data->elastic_fusion_xforms );
    mean[3] = calculate_mean( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              seq_data->elastic_fusion_xforms );
    sdev[3] = calculate_stddev( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                                mean[3], seq_data->elastic_fusion_xforms );
  }

  rmse[4] = NAN; mean[4] = NAN; sdev[4] = NAN;
  if ( seq_data->kintinuous_xforms )
  {
    rmse[4] = calculate_rmse( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              seq_data->kintinuous_xforms );
    mean[4] = calculate_mean( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                              seq_data->kintinuous_xforms );
    sdev[4] = calculate_stddev( seq_data->gt_corrs, seq_data->n_gt_corrs, 
                                mean[4], seq_data->kintinuous_xforms );
  }
}
}

////////////////////////////////////////////////////////////////////////////////
// MiscMath
////////////////////////////////////////////////////////////////////////////////
namespace mmath
{
template< typename T>
static inline T
linmap( T val, T a_1, T a_2, T b_1 = 0.0, T b_2 = 1.0 )
{
  if ( val < a_1 ) return b_1;
  if ( val > a_2 ) return b_2;
  T a = ( b_1 - b_2 ) / ( a_1 - a_2 );
  T b = b_2 - a_2 * a;
  return (a * val + b);
}

template <typename T>
static inline std::vector<T>
linspace( T a, T b, int n_elem = 100 )
{
  assert( a != b );
  T diff = b - a;
  T range = std::abs( diff );
  int sign = static_cast<int>( diff / range );
  double step = static_cast<double>(range) / 
                static_cast<double>( n_elem - 1 ) * sign;
  std::vector<T> linspace;
  linspace.reserve( n_elem );
  for ( int i = 0 ; i < n_elem ; ++i )
  {
    linspace.push_back( a + i * step );
  }
  return linspace;
}

template <typename T>
static inline T
gauss( T x, T mu, T sigma )
{
    return expf( - ( ( x - mu ) * ( x - mu ) ) / ( (T)2 * sigma * sigma ) );
}

template <typename T>
static inline std::vector<T>
gaussmf( std::vector<T> vals, T mu, T sigma )
{
  T two_sigma_sq = (T)2 * sigma * sigma;
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    vals[i] = std::exp( - ( ( vals[i] - mu ) * ( vals[i] - mu ) ) / two_sigma_sq );
  }
  return vals;
}

static inline float
sigmoid( float x, float mean, float sigma )
{
  return 1.0f / (1.0f + expf( -(sigma * ( x - mean ) ) ) );
}

static inline std::vector<float>
pow_vector( std::vector<float> vals, float k )
{
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    vals[i] = pow( vals[i], k );
  }
  return vals;
}

static inline std::vector<float>
exp_vector( std::vector<float> vals, float k )
{
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    vals[i] = exp( - k * vals[i] );
  }
  return vals;
}

static inline std::vector<float>
sin_vector( std::vector<float> vals)
{
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    vals[i] = sin( vals[i] );
  }
  return vals;
}

static inline std::vector<float>
cos_vector( std::vector<float> vals )
{
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    vals[i] = cos( vals[i] );
  }
  return vals;
}

static float
clamp( r32 cur, r32 init, r32 fin )
{
  r32 res;
  if( init < fin )
  {
    res = bsc::min( cur, fin );
  }
  else
  {
    res = bsc::max( cur, fin );
  }
  return res;
}

}

template <typename T>
static inline std::vector<T>
operator+( std::vector<T> vals, T val )
{
  std::vector<T> results;
  results.reserve( vals.size() );
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    results.push_back( vals[i] + val );
  }
  return results;
}

// NOTE: Operator overloads need to be in the global namespace for them to work.
// Dislike this.
template <typename T>
static inline std::vector<T>
operator+( T val, std::vector<T> vals )
{
  std::vector<T> results;
  results.reserve( vals.size() );
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    results.push_back( vals[i] + val );
  }
  return results;
}

template <typename T>
static inline std::vector<T>
operator*( T val, std::vector<T> vals )
{
  std::vector<T> results;
  results.reserve( vals.size() );
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    results.push_back( vals[i] * val );
  }
  return results;
}

template <typename T>
static inline std::vector<T>
operator*( std::vector<T> vals, T val )
{
  std::vector<T> results;
  results.reserve( vals.size() );
  for ( int i = 0 ; i < vals.size() ; ++i )
  {
    results.push_back( vals[i] * val );
  }
  return results;
}


////////////////////////////////////////////////////////////////////////////////
// Helper functions with no better place to be (?)
////////////////////////////////////////////////////////////////////////////////

void 
filter_reconstruction( FETReconstruction * reconstruction,
                      const OptimizationParameters * optim_params,
                      bool remove_faraway = false )
{
  std::vector<FETFeature*> features_to_remove;
  for ( i32 shape_idx = 0; shape_idx < reconstruction->NShapes(); shape_idx++ )
  {
    FETShape * shape = reconstruction->Shape( shape_idx );
    for ( i32 feat_idx = 0; feat_idx < shape->NFeatures(); feat_idx++ )
    {
      FETFeature * f = shape->Feature( feat_idx );

      // make salience linear function of depth, simply;
      double depth = -f->Position().Z();
      double salience = mmath::linmap( depth, 1.0, optim_params->max_depth + 0.01, 1.0, 0.0 );
      f->SetSalience( salience );
    }

    // Create kdtrees ahead of time
    if (!shape->kdtree) 
    {
      FETFeature tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
      shape->kdtree = new R3Kdtree<FETFeature *>(shape->features, position_offset);
      if (!shape->kdtree) RNAbort("Cannot build kdtree");
    }

  }

  for ( i32 i = 0; i < features_to_remove.size() ; ++i )
  {
    delete features_to_remove[i];
  }
}

void
calculate_parametrization( SequenceData * seq_data,
                           const Options * opts )
{
  std::vector<double> distances( opts->frame_end, 0.0 );
  seq_data->parametrization.resize( opts->frame_end, 0.0 );
  
  // Create gaussian kernel
  int kernel_size = 7;
  std::vector<double> vals = mmath::linspace( 0.0, (double)kernel_size-1, kernel_size );
  std::vector<double> gauss_kernel = mmath::gaussmf( vals, (kernel_size-1.0)/2.0, 1.5 );

  /* compute parametrization along the path */
  for( int i = 1; i < opts->frame_end; i++ )
  {
    bsc::vec3d origin0 = bsc::vec3d(seq_data->current_xforms[i-1][3]);
    bsc::vec3d origin1 = bsc::vec3d(seq_data->current_xforms[i][3]);
    double d = bsc::norm(origin0 - origin1);
    distances[i] = distances[i - 1] + d;
  }

  /* smooth out the parametrization values */
  for( int frame_idx = 0; frame_idx < opts->frame_end; frame_idx++ )
  {
    seq_data->parametrization[frame_idx] = 0.0;
    double sum = 0.0;
    for( int kernel_idx = -kernel_size/2 ; kernel_idx <= kernel_size/2 ; kernel_idx++ )
    {
      int idx = frame_idx + kernel_idx;
      if ( idx >= 0 && idx < opts->frame_end ) 
      {
        double kernel_val = gauss_kernel[kernel_idx+(kernel_size/2)]; 
        sum += kernel_val;
        seq_data->parametrization[frame_idx] += distances[idx] * kernel_val;
      }
    }
    seq_data->parametrization[frame_idx] /= sum;
  }
}

void
update_bboxes( SequenceData * seq_data, const Options * opts, float max_dist )
{
  FETReconstruction *r = seq_data->reconstruction;
  if ( !seq_data->bboxes ) 
  {
    seq_data->bboxes = new R3Box[ opts->frame_end ]; 
  }

  for( i32 i = 0; i < opts->frame_end; i++ )
  {
    FETShape *s = r->Shape(i);
    R3Affine t( bscmat4_to_R4Matrix(seq_data->current_xforms[i] ) );
    seq_data->bboxes[i] = R3null_box;
    for( int j = 0 ; j < s->NFeatures() ; ++j )
    {
      R3Point p = s->Feature(j)->Position();
      if( -p.Z() < max_dist )
      {
        p.Transform(t);
        seq_data->bboxes[i].Union( p );
      }
    }
  }
}

void 
communicate_changes( const FETReconstruction * reconstruction,
                     SequenceData * seq_data, 
                     OptimizationParameters * optim_params,
                     const Options * opts )
{
  for (int i = 0; i < opts->frame_end; i += 1)
  {
    bsc::mat4d xform = R4Matrix_to_bscmat4(
        reconstruction->Shape(i)->current_transformation.Matrix());
    seq_data->xforms_history[seq_data->n_xforms * seq_data->n_frames + i] =
        xform;
  }

  // Update current transformations pointer
  seq_data->current_xforms = &(seq_data->xforms_history[seq_data->n_xforms * seq_data->n_frames]);

  seq_data->n_xforms += 1;
  optim_params->iter_idx += 1;
}

void
update_shape_proxy_map( StructuralModel * model )
{
  for ( i32 i = 0 ; 
        i < model->structure->NProxies() ; 
        ++i )
  {
    
    Proxy * proxy = model->structure->GetProxy( i);
    model->shape_to_proxy->insert( { proxy->shape, proxy } );
  }
}

void write_corr_vs_dist( const char * filename, 
                         const Correspondence * corrs,
                         const i32 n_corrs,
                         const bsc::mat4d * xforms )
{
  FILE * f = fopen( filename, "w" );
  if ( !f ) 
  {
    printf("Failed to open the file %s!\n", filename );
    return;
  }

  for ( i32 corr_idx = 0;
        corr_idx < n_corrs;
        ++corr_idx )
  {
    const Correspondence * corr = &(corrs[corr_idx]);

    bsc::vec4d pt_A = bsc::vec4d( corr->point_A, 1.0 );
    bsc::vec4d pt_B = bsc::vec4d( corr->point_B, 1.0 );

    pt_A = xforms[corr->idx_A] * pt_A;
    pt_B = xforms[corr->idx_B] * pt_B;

    r32 cur_err  = bsc::norm( bsc::vec3d(pt_A) - bsc::vec3d(pt_B) );
    r32 cur_dist =  abs( corr->idx_A - corr->idx_B );
    fprintf(f, "%f %f\n", cur_err, cur_dist );
  }
  
  fclose(f);
}

////////////////////////////////////////////////////////////////////////////////
// Rigid alignment
////////////////////////////////////////////////////////////////////////////////
namespace align
{
void
reconstruction_to_world ( SequenceData * seq_data, 
                          const Options *opts,
                          const i32 n_iter = 1,
                          i32 n_frames = -1,
                          bsc::mat4d * xforms_ptr = NULL )
{
  bsc::stop_watch timer;
  timer.start();

  FETReconstruction * recstr = seq_data->reconstruction;
  bool update_shapes = false;
  if ( xforms_ptr == NULL )
  {
    update_shapes = true;
    xforms_ptr = seq_data->current_xforms;
  }
  if( n_frames == -1 )
  {
    n_frames = opts->frame_end;
  }

  bsc::mat4d align;

  for ( i32 iter = 0 ; iter < n_iter ; ++iter )
  {
    
    // find inliers of vectors in each dimension + centroid
    R3Vector x = R3zero_vector;
    R3Vector y = R3zero_vector;
    R3Vector z = R3zero_vector;
    R3Point centroid = R3zero_point;
    r32 x_count = 0, y_count = 0, z_count = 0, c_count = 0;
    for ( i32 shape_idx = 0 ; shape_idx < n_frames ; shape_idx += 1 )
    {
      FETShape * shape = recstr->Shape( shape_idx );
      R3Affine t( bscmat4_to_R4Matrix( xforms_ptr[shape_idx ]) );

      for ( i32 feat_idx = 0 ;
            feat_idx < shape->NFeatures() ;
            feat_idx++ )
      {
        FETFeature * feat = shape->Feature( feat_idx );
        R3Point pos = feat->Position();

        bool is_nan = std::isnan( pos.X() ) ||
                      std::isnan( pos.Y() ) ||
                      std::isnan( pos.Z() );
        if( is_nan )
        {
          printf("Shape_idx %d Feat_idx %d %f %f %f\n",
          shape_idx, feat_idx, pos.X(), pos.Y(), pos.Z() );
        }

        if ( pos.Z() < -2.5 || is_nan ) continue;
        
        pos.Transform(t);

        centroid += pos;
        c_count++;

        R3Vector normal = feat->Normal();
        normal.Transform(t);
        normal.Normalize();
        r32 angle_x = std::min( fabs(R3InteriorAngle(normal, R3posx_vector)),
                                fabs(R3InteriorAngle(normal, R3negx_vector)) );
        r32 angle_y = std::min( fabs(R3InteriorAngle(normal, R3posy_vector)),
                                fabs(R3InteriorAngle(normal, R3negy_vector)) );
        r32 angle_z = std::min( fabs(R3InteriorAngle(normal, R3posz_vector)),
                                fabs(R3InteriorAngle(normal, R3negz_vector)) );
        r32 min_angle = bsc::min3( angle_x, angle_y, angle_z );
        if ( min_angle == angle_x )
        {
          if ( normal.Dot( R3posx_vector) < 0 ) normal = -normal;
          x += normal;
          x_count += 1; 
        }
        if ( min_angle == angle_y )
        {
          if ( normal.Dot( R3posy_vector) < 0 ) normal = -normal;
          y += normal;
          y_count += 1; 
        }
        if ( min_angle == angle_z )
        {
          if ( normal.Dot( R3posz_vector) < 0 ) normal = -normal;
          z += normal;
          z_count += 1; 
        }
      
      }
    }

    x /= x_count; x.Normalize();
    y /= y_count; y.Normalize();
    z /= z_count; z.Normalize();
    centroid /= c_count;

    // find vector with most inliers -> we will map it to up vector
    r32 counts[3]  = { x_count, y_count, z_count };
    R3Vector vecs[3] = { x, y, z };
    i32 ind[3] = { 0, 1, 2 };

    // sorting with indices trick
    // TODO(maciej): Generalize this using c++11 maybe?
    array_helper = &(counts[0]);
    qsort( ind, 3, sizeof(ind[0]), arr_idx_compare );

    // create orthonormal basis
    vecs[ ind[0] ] = vecs[ ind[2] ];
    vecs[ ind[0] ].Cross( vecs[ ind[1] ]);
    vecs[ ind[0] ].Normalize();

    vecs[ ind[1] ] = vecs[ ind[2] ];
    vecs[ ind[1] ].Cross( vecs[ ind[0] ] );
    vecs[ ind[1] ].Normalize();

    // construct the matrix where most prominent axis maps to Y(up vector)
    align[0] = bsc::vec4d( vecs[ ind[0] ].X(), 
                           vecs[ ind[0] ].Y(), 
                           vecs[ ind[0] ].Z(), 
                           0.0 );
    align[1] = bsc::vec4d( vecs[ ind[2] ].X(), 
                           vecs[ ind[2] ].Y(), 
                           vecs[ ind[2] ].Z(), 
                           0.0 );
    align[2] = bsc::vec4d( vecs[ ind[1] ].X(), 
                           vecs[ ind[1] ].Y(), 
                           vecs[ ind[1] ].Z(), 
                           0.0 );
    align[3] = bsc::vec4d( centroid.X(),
                           centroid.Y(),
                           centroid.Z(),
                           1.0 );
    align = bsc::inverse(align);


    // apply transformations
    for (i32 i = 0;  i < n_frames; i++ )
    {
      xforms_ptr[i] = align * xforms_ptr[i];
    }

    if( update_shapes ) 
    {
      for( i32 i = 0; i < recstr->NShapes() ; i++ )
      {
        FETShape * shape = recstr->Shape( i );
        bsc::mat4d xform;
        if( i < n_frames ) { xform = xforms_ptr[i]; }
        else 
        {
          bsc::mat4d m = R4Matrix_to_bscmat4(shape->current_transformation.Matrix());
          xform = align * m;
        }
        shape->current_transformation = R3Affine(bscmat4_to_R4Matrix(xform));
      }
    }


  }

  // compute average camera basis
  bsc::vec3d avg_n( 0.0, 0.0, 0.0 );
  bsc::vec3d avg_u( 0.0, 0.0, 0.0 );

  for (i32 i = 0; i < n_frames;  i++ )
  {
    avg_n += bsc::vec3d( xforms_ptr[i][1]);
    avg_u += bsc::vec3d( xforms_ptr[i][0]);
  }
  avg_n /= n_frames;
  avg_u /= n_frames;
  avg_n = bsc::normalize( avg_n );
  avg_u = bsc::normalize( avg_u );

  // figure out what is best approaximating basis based on 
  // highest compnent in each average
  r64 max_u = bsc::max3( avg_u.x, avg_u.y, avg_u.z );
  r64 min_u = bsc::min3( avg_u.x, avg_u.y, avg_u.z );
  r64 top_u = fmax( fabs(min_u), max_u );
  r64 sign_u = (fabs(min_u) > max_u) ? -1.0 : 1.0;
  top_u *= sign_u;
  bsc::vec3d u( (top_u == avg_u.x) * sign_u, 
              (top_u == avg_u.y) * sign_u, 
              (top_u == avg_u.z) * sign_u);
  
  r64 max_n = bsc::max3( avg_n.x, avg_n.y, avg_n.z );
  r64 min_n = bsc::min3( avg_n.x, avg_n.y, avg_n.z );
  r64 top_n = fmax( fabs(min_n), max_n );
  r64 sign_n = (fabs(min_n) > max_n) ? -1.0f : 1.0f;
  top_n *= sign_n;
  bsc::vec3d n( (top_n == avg_n.x) * sign_n, 
              (top_n == avg_n.y) * sign_n, 
              (top_n == avg_n.z) * sign_n);
  
  if ( fabs(bsc::dot(n, u)) - 1.00 <= 0.0001f )
  {
    u = bsc::vec3d( (top_n != avg_n.x) * -sign_n, 
                  (top_n != avg_n.y) * -sign_n, 
                  (top_n != avg_n.z) * -sign_n);
    for ( int i = 0 ; i < 3 ; ++i )
    {
      if ( fabs(u[i]) == 1.0 )
      {
        u[i] = 0.0f;
        break;
      }
    }
  }

  bsc::mat4d post_align;
  post_align[0] = bsc::vec4d(u, 0.0f);
  post_align[1] = bsc::vec4d(n, 0.0f);
  post_align[2] = bsc::vec4d( bsc::normalize( bsc::cross( u, n ) ), 0.0f);

  // fix current
  bsc::mat4d inv_post_align = bsc::inverse(post_align);
  for (i32 i = 0; i < n_frames; i++)
  {
    bsc::mat4d xform = inv_post_align * xforms_ptr[i];
    xforms_ptr[i] = xform;
  }

  if( update_shapes ) 
  {
    for( i32 i = 0; i < recstr->NShapes() ; i++ )
    {
      FETShape * shape = recstr->Shape( i );
      bsc::mat4d xform;
      if( i < n_frames ) { xform = xforms_ptr[i]; }
      else 
      {
        bsc::mat4d m = R4Matrix_to_bscmat4(shape->current_transformation.Matrix());
        xform = inv_post_align * m;
      }
      shape->current_transformation = R3Affine(bscmat4_to_R4Matrix(xform));
    }
  }
  if( opts->print_verbose ) printf("Alignment with world took %f sec.\n", timer.elapsed() );
}
}
