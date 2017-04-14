// FineToCoarseRegistration - fetregister
//
// This is the main file of the project - it reads in the required files and 
// sets up the optimization. If you wish to modify this code base you should 
// most likely start with initialize_optimization(...) and 
// optimize_camera_poses(...) functions.

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

#define FETREGISTER

#include "common.h"
#define IO_IMPLEMENTATION
#include "io.h"
#define OPTIMIZATION_IMPLEMENTATION
#include "optimization.h"

#ifdef FETREGISTER_USE_WINDOW
#define RENDERING_IMPLEMENTATION
#include "rendering.h"
#define GUI_IMPLEMENTATION
#include "gui.h"
#endif


////////////////////////////////////////////////////////////////////////////////
// Global State and Forward Decl
////////////////////////////////////////////////////////////////////////////////
static Options opts;
static SequenceData seq_data;
static ViewOptions view_opts;
static OptimizationParameters optim_params;

#ifdef FETREGISTER_USE_WINDOW
static RenderingState render;
#endif
////////////////////////////////////////////////////////////////////////////////
// Optimization Function
////////////////////////////////////////////////////////////////////////////////

void 
optimize_camera_poses( SequenceData *seq_data,
                       OptimizationParameters *optim_params,
                       const Options *opts)
{
  RNTime t, t_all;
  t_all.Read();

  // Info 
  if( opts->print_verbose )
  {
    printf("\nIteration %d\n", optim_params->iter_idx);
    printf("   Weights    -> CP: %f ; Trajectory : %f ; "
           "Inertia : %f ; Geometry: %f ; Hierarchy : %f\n",
            optim_params->cp_weight,
            optim_params->traj_weight,
            optim_params->inertia_weight,
            optim_params->geometrical_weight,
            optim_params->hierarchical_weight );
  }

  // Initial setup and cleanup 
  FETReconstruction *recstr = seq_data->reconstruction;

  optcorrs::clear_correspondences( seq_data->cp_corrs );
  optcorrs::clear_correspondences( seq_data->geometrical_corrs );
  optcorrs::clear_correspondences( seq_data->hierarchical_corrs );

  // we want to increment exponent for the decreasing weights only once per iteration
  optproxy::reset_seen_index_increments( &(optim_params->seen_count) );
  recstr->discard_outliers       = true;
  recstr->discard_boundaries     = true;
  if( optim_params->geometrical_weight > 0 && 
      optim_params->hierarchical_weight > 0 )
  {
   optstruct::create_hierarchical_structure( 0, opts->frame_end, 
                                                 seq_data, optim_params, opts );
  }

  // Control parametrization increase
  static int renormalized = false;
  if ( optim_params->segment_length >= 2 * seq_data->parametrization[opts->frame_end-1] &&
       !renormalized )
  {
    renormalized = true;
    optim_params->segment_length = 2 * seq_data->parametrization[opts->frame_end-1];
  } 
  // Realistically we should not be longer than 100km
  if ( optim_params->segment_length >= 100000 ) 
  {
    optim_params->segment_length = 100000;
  }

  // closest point correspondence insertion
  if( optim_params->cp_weight > 0.0f)
  {
    optcorrs::add_closest_point_correspondeces( recstr, seq_data, 
                                                optim_params, opts );
  }

  // Structure creation
  if( optim_params->geometrical_weight > 0 && 
      optim_params->hierarchical_weight > 0 )
  {
    if( opts->print_verbose ) printf("\n  Geometrical Relationships: \n" );
    optstruct::create_hierarchical_structure( 1, opts->frame_end,
                                              seq_data,
                                              optim_params, 
                                              opts );
    if( opts->print_verbose )  printf("      -- Created %d proxies\n", 
                                       seq_data->model->structure->NProxies() );
  }


  /* Structure is done -> Add correspondences! */
  if( optim_params->geometrical_weight > 0.0f )
  {
   optstruct::create_geometrical_relationships( 1, seq_data, optim_params );
   if( opts->print_verbose ) printf("      -- Geometrical Relationships : %lu\n", 
                                       seq_data->model->relationships->size() );
   optstruct::add_geometrical_correspondeces( seq_data->model->relationships, 
                                              seq_data->geometrical_corrs,
                                              opts );
  }

  if( optim_params->hierarchical_weight > 0.0f )
  {
    optstruct::add_hierarchical_correspondences( seq_data->reconstruction,
                                                seq_data->model,
                                                seq_data->hierarchical_corrs, 
                                                opts );
  }

  // Optimization call
  t.Read();
  if( opts->print_verbose ) printf("\n  Starting Optimization for "
                                   "iteratrion: %d:\n", optim_params->iter_idx);


  optimize_transformations(recstr, seq_data, optim_params, opts);

  if( opts->print_verbose ) printf("  done in %f sec.\n", t.Elapsed());

  
  // COMMUNICATION -> Infroms the rest of the world on what happened
  communicate_changes(recstr, seq_data, optim_params, opts );
  align::reconstruction_to_world( seq_data, opts, 4);
  align::reconstruction_to_world( seq_data, opts, 4);
  update_bboxes( seq_data, opts, optim_params->max_depth );
  calculate_parametrization( seq_data, opts );
  view_opts.xform_idx  = seq_data->n_xforms - 1;
  view_opts.xforms_ptr = &(seq_data->xforms_history[ 
                                   view_opts.xform_idx * seq_data->n_frames ] );
  seq_data->errors[seq_data->n_xforms - 1] =
                     eval::calculate_rmse(seq_data->gt_corrs,
                                          seq_data->n_gt_corrs,
                                          seq_data->current_xforms,
                                          seq_data->per_correspondence_error );
  // Report status
  if( opts->print_verbose )
  {
    printf("\nIteration %d done in %f sec. Current error %f \n", 
                                        optim_params->iter_idx - 1, 
                                        t_all.Elapsed(),
                                        seq_data->errors[seq_data->n_xforms - 1]);
  }
}


i32 initialize_optimization( SequenceData *seq_data,
                            OptimizationParameters *optim_params,
                            Options *opts )
{
  if( !io::read_configuration_file( opts->conf_filename,
                                    seq_data,
                                    opts ) )
  {
    return 0;
  }

  if( !io::read_fet_file( opts->fet_filename,
                          seq_data,
                          opts ) )
  {
    return 0;
  }


  // apply the pairwise transformatios to our reconstruction
  if( seq_data->pairwise_xforms != NULL )
  {
    for (i32 i = 0; i < opts->frame_end; i++ )
    {
      seq_data->current_xforms[i] = seq_data->pairwise_xforms[i];
      seq_data->initial_xforms[i] = seq_data->pairwise_xforms[i];
    }
  }
  else
  {
    printf( "Please provide pairwise matches in the configuration file.\n" );
  }

  // remove possibly invalid corrs ( issue with RR )
  eval::invalidate_corrs( seq_data, seq_data->robust_reconstruction_xforms );

  // pass through features in fet file, remove unwanted ones
  filter_reconstruction( seq_data->reconstruction, optim_params );
  seq_data->reconstruction->solver = opts->solver;

  // align the transformations and copy them relevant storage
  align::reconstruction_to_world( seq_data, opts, 4);
  align::reconstruction_to_world( seq_data, opts, 4);

  if( seq_data->robust_reconstruction_xforms != NULL)
  {
    // the robust reconstruction matrices are not aligned with world
    align::reconstruction_to_world( seq_data, opts, 4, opts->frame_end, 
                                    seq_data->robust_reconstruction_xforms );
  }

  if( seq_data->elastic_fusion_xforms != NULL)
  {
    // the elastic fusion matrices are not aligned with world
    align::reconstruction_to_world( seq_data, opts, 4, opts->frame_end, 
                          seq_data->elastic_fusion_xforms );
  }

  if( seq_data->kintinuous_xforms != NULL)
  {
    // the kintinuous matrices are not aligned with world
    align::reconstruction_to_world( seq_data, opts, 4, opts->frame_end, 
                          seq_data->kintinuous_xforms );
  }

  seq_data->errors[seq_data->n_xforms - 1] =
                     eval::calculate_rmse(seq_data->gt_corrs,
                                          seq_data->n_gt_corrs,
                                          seq_data->current_xforms,
                                          seq_data->per_correspondence_error );
  update_bboxes( seq_data, opts, optim_params->max_depth );
  calculate_parametrization( seq_data, opts );

  optim_params->cp_weight           = optim_params->init_cp_weight;      
  optim_params->traj_weight         = optim_params->init_traj_weight;    
  optim_params->geometrical_weight  = optim_params->init_geometrical_weight;    
  optim_params->hierarchical_weight = optim_params->init_hierarchical_weight;  
  optim_params->inertia_weight      = optim_params->init_inertia_weight;      
  optim_params->angle_threshold     = optim_params->init_angle_threshold;
  optim_params->distance_threshold  = optim_params->init_distance_threshold;
  optim_params->max_n_corrs         = opts->frame_end * 50;

  return 1;
}

void
run_optimization( SequenceData * seq_data,
                 OptimizationParameters * optim_params,
                 ViewOptions * view_opts,
                 const Options * opts )
{
  RNTime t1;
  t1.Read();

  // Weights
  r32 cp_weight_step       = (optim_params->final_cp_weight - 
                              optim_params->init_cp_weight) / 
                              (optim_params->n_iter - 1);
  r32 traj_weight_step     = (optim_params->final_traj_weight - 
                              optim_params->init_traj_weight) / 
                              (optim_params->n_iter - 1);
  r32 inertia_weight_step  = (optim_params->final_inertia_weight - 
                              optim_params->init_inertia_weight) / 
                              (optim_params->n_iter - 1);
  r32 geometrical_weight_step   = (optim_params->final_geometrical_weight - 
                              optim_params->init_geometrical_weight) / 
                              (optim_params->n_iter - 1);
  r32 hierarchical_weight_step = (optim_params->final_hierarchical_weight - 
                              optim_params->init_hierarchical_weight) / 
                              (optim_params->n_iter - 1);

  for ( i32 idx = 0 ; 
        idx < optim_params->n_iter;
        idx++ )
  {
    optimize_camera_poses( seq_data, optim_params, opts );
    optim_params->segment_length *= optim_params->segment_length_growth;

    if( optim_params->n_iter > 1 )
    {
      optim_params->cp_weight       += cp_weight_step;
      optim_params->traj_weight     += traj_weight_step;
      optim_params->inertia_weight  += inertia_weight_step;
      optim_params->hierarchical_weight += hierarchical_weight_step;
      optim_params->geometrical_weight   += geometrical_weight_step;

      optim_params->cp_weight = mmath::clamp(optim_params->cp_weight,
                                                 optim_params->init_cp_weight,
                                                 optim_params->final_cp_weight);
      optim_params->traj_weight = mmath::clamp(optim_params->traj_weight,
                                               optim_params->init_traj_weight,
                                               optim_params->final_traj_weight);
      optim_params->inertia_weight = mmath::clamp(optim_params->inertia_weight,
                                            optim_params->init_inertia_weight,
                                            optim_params->final_inertia_weight);
      optim_params->hierarchical_weight = mmath::clamp(optim_params->hierarchical_weight,
                                           optim_params->init_hierarchical_weight,
                                           optim_params->final_hierarchical_weight);
      optim_params->geometrical_weight = mmath::clamp(optim_params->geometrical_weight,
                                             optim_params->init_geometrical_weight,
                                             optim_params->final_geometrical_weight);
    }
  }

  // make sure we are showing most current transformation after optimization
  view_opts->xform_idx     = seq_data->n_xforms - 1;
  view_opts->xforms_ptr    = &(seq_data->xforms_history[ 
                               view_opts->xform_idx * seq_data->n_frames ] );

  // write errors
  if( opts->output_error_filename )
  {
    io::write_errors( opts->output_error_filename, seq_data );
  }
  if( opts->print_verbose )
  {
    printf("Optimization done in %f sec. Current error: %f \n", 
                       t1.Elapsed(), seq_data->errors[seq_data->n_xforms - 1] );
  }
}

////////////////////////////////////////////////////////////////////////////////
// Shortcuts
////////////////////////////////////////////////////////////////////////////////
#ifdef FETREGISTER_USE_WINDOW

void keyboard(ViewOptions *view_opts,
              RenderingState *render,
              const SequenceData *seq_data,
              const Options *opts)
{
  if( bsc::ui::is_key_pressed(bsc::ui::key_r))
  {
    if( view_opts->xform_idx < seq_data->n_xforms - 1)
    {
      view_opts->xforms_ptr = seq_data->current_xforms;
      view_opts->xform_idx = seq_data->n_xforms - 1;
    }
    else
    {
      view_opts->xforms_ptr = seq_data->initial_xforms;
      view_opts->xform_idx = 0;
    }

    view_opts->xform_scheme = (view_opts->xform_scheme == CURRENT_TRANSFORMATION) ? 
                               INITIAL_TRANSFORMATION:
                               CURRENT_TRANSFORMATION;

    rendering::request_update( render, opts );
  }
  if( bsc::ui::is_key_pressed(bsc::ui::key_i))
  {
    view_opts->inspection_mode = !view_opts->inspection_mode;
  }
  if( bsc::ui::is_key_pressed(bsc::ui::key_b))
  {
    view_opts->show_bboxes = !view_opts->show_bboxes;
  }

  if( bsc::ui::is_key_pressed(bsc::ui::key_c))
  {
    view_opts->show_cp_corrs = !view_opts->show_cp_corrs;
  }

  if( bsc::ui::is_key_pressed(bsc::ui::key_d))
  {
    view_opts->feat_disp_mode = (FeatureDisplayMode)((view_opts->feat_disp_mode + 1) % 5);
  }

  if( bsc::ui::is_key_pressed(bsc::ui::key_t))
  {
    view_opts->show_trajectory = !view_opts->show_trajectory;
  }

  if( bsc::ui::is_key_pressed(bsc::ui::key_s))
  {
    view_opts->show_frames = !view_opts->show_frames;
  }

  if( opts->output_conf_filename )
  {
    if( bsc::ui::is_key_pressed(bsc::ui::key_w))
    {
      io::write_configuration_file( opts->output_conf_filename,
                                    seq_data,
                                    opts);
      if( strcmp( seq_data->dataset, "tum" ) == 0 )
      {
        io::write_tum_trajectory_file( "tum.traj", seq_data, opts );
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Init and Draw Functions
////////////////////////////////////////////////////////////////////////////////

static i64 total_n_points = 0;
i32 Init()
{
  bsc::stop_watch timer;
  timer.start();

  if( !opts.load_frames )
  {
    if( !initialize_optimization( &seq_data, &optim_params, &opts ) )
    {
      return 0;
    }
  }
  else
  {
    if( !io::read_configuration_file( opts.conf_filename,
                                      &seq_data,
                                      &opts))
    {
      return 0;
    }
    if( seq_data.pairwise_xforms != NULL )
    {
      for (i32 i = 0; i < opts.frame_end; i++ )
      {
        seq_data.current_xforms[i] = seq_data.pairwise_xforms[i];
        seq_data.initial_xforms[i] = seq_data.pairwise_xforms[i];
      }
    }
  }

  // Viewing Initialization  
  // Only available after parsing arguments
  view_opts.frame_idx_B = opts.frame_end - 1;
  view_opts.xforms_ptr = seq_data.initial_xforms;

   // set some options only available after reading the conf
  render.view = render.cam_controls.initialize(&render.cam,
                                               bsc::vec3(5.0f, 5.0f, 5.0f),
                                               bsc::vec3(0.0f, 0.0f, 0.0f),
                                               bsc::vec3(0.0f, 1.0f, 0.0f));

  // initialize rendering -> common shaders and geometries
  bsc::vec4i viewport(0, 0, bsc::g_window.size.x - 256, bsc::g_window.size.y);
  rendering::init(&render,
                  &seq_data,
                  &opts,
                  viewport);

  // initialize pointcloud rendering ->
  if( opts.load_frames )
  {
    total_n_points = rendering::init_pointclouds_from_frames(&render,
                                                             &seq_data,
                                                             &opts);
  }
  else
  {
    total_n_points = rendering::init_pointclouds_from_reconstruction(&render,
                                                                     &seq_data,
                                                                     &view_opts,
                                                                     &opts);
  }
  if( opts.print_verbose)
  {
    printf("-----------------------------\n");
    printf("Initialization took %f sec.\n", timer.elapsed());
  }

  return 1;
}

void Display()
{
  glDisable(GL_CULL_FACE);

  // Update IMGUI and user interface
  bsc::vec4i viewport(0, 0, bsc::g_window.size.x - 256, bsc::g_window.size.y);
  bsc::ui::begin_frame();

  keyboard(&view_opts, &render, &seq_data, &opts);

  if( bsc::ui::is_mouse_over(viewport) && !ImGui::IsMouseHoveringAnyWindow())
  {
    render.cam_controls.update(&render.cam, &render.view, &viewport);
  }

  //TODO: (maciej)Write control wrapper for the projection matrices
  static r32 near = 0.1f;
  static r32 far = 100.0f;
  r32 fovy = bsc::deg2rad(view_opts.fov_y);
  r32 aspect_ratio = (r32)viewport.z / viewport.w;
  render.proj = bsc::perspective(fovy, aspect_ratio, near, far);

  // Main drawing call
  rendering::draw(&render, &seq_data, &view_opts, 
                  &optim_params, &opts, viewport);
  
  // Inspection code
  // TODO: Move it to relative place. GUI?
  if( view_opts.inspection_mode )
  {
    i32 lmb_clicked = bsc::ui::is_mouse_clicked( bsc::ui::lmb, bsc::ui::ctrl ); 
    
    // Read values
    bsc::ui::state & state = bsc::get_UI();
    bsc::vec2d pos = bsc::vec2d( bsc::g_window.pixel_ratio * state.mouse_pos.x,
                            bsc::g_window.pixel_ratio * (bsc::g_window.size.y - state.mouse_pos.y) );
    float color[4] = {-1.0, -1.0, -1.0, -1.0};
    glBindFramebuffer(GL_READ_FRAMEBUFFER, render.fullscreen_framebuffer);
    glReadBuffer( GL_COLOR_ATTACHMENT1 );
    glReadPixels(	pos.x, pos.y, 1, 1, 
                  GL_RG, GL_FLOAT, (GLvoid*)&(color[0]) );
    glReadBuffer(GL_NONE);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    int class_id = color[0];
    int object_id = color[1];
    
    // Chose action based on read out values
    if( class_id > 0 )
    {
      if( lmb_clicked )
      {
        view_opts.selected_class_id = class_id;
        view_opts.selected_object_id = object_id;
        render.geometrical_corrs_dirty = true;
      }

      ImGui::SetNextWindowPos(ImVec2(state.mouse_pos.x+10, state.mouse_pos.y+10));
      ImGui::SetNextWindowSize(ImVec2(264, 164));
      ImGui::Begin("Info", NULL, ImVec2(-1, -1), -1.0f,
                ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                    ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
      if( class_id == POINTCLOUD_ID )
      {
          ImGui::Text("Pointcloud");
          ImGui::Text(" - Frame %d", object_id);
          ImGui::Text(" - SIFT: %d", seq_data.pairwise_xforms_inliers[object_id] );
      }
      else if( class_id == PROXY_ID )
      {
        ImGui::Text("Proxy");
        Proxy * p = seq_data.model->structure->GetProxy(object_id);
        ImGui::Text("-ShapeIDX %d %d", p->shape->reconstruction_index, p->shape->NFeatures() );
        ImGui::Text("-IDX %d %d", object_id, p->structure_index );
        ImGui::Text("-Level %d", p->level );
        double dist = seq_data.parametrization[p->max_shape_index] -
                      seq_data.parametrization[p->min_shape_index];
        ImGui::Text("-Limits %d %d | %f", p->min_shape_index, p->max_shape_index, dist );
        ImGui::Text("-Inliers %d", p->n_inliers );
        ImGui::Text("-Salience %f", p->feature->salience );
        ImGui::Text("-Stats: Mean %f SDev %f", p->mean, p->sdev );
        R3Vector normal = p->feature->Normal(TRUE);
        ImGui::Text("-Normal %f %f %f\n", normal.X(), normal.Y(), normal.Z());
      }
      else if( class_id == CP_CORR_ID )
      {
          ImGui::Text("Cp_corr");
      }
      else if( class_id == HIERARCHICAL_CORR_ID )
      {
          ImGui::Text("hierarchical_corr");
      }
      else if( class_id == GEOMETRICAL_CORR_ID )
      {
        ImGui::Text("Structural Corr.");
        ImGui::Text("Object ID: %d/%lu\n", object_id, seq_data.geometrical_corrs.size());
        // Unpack some data
        FETCorrespondence * corr = seq_data.geometrical_corrs[ object_id ];
        if ( corr->RelationshipType() == PARALLEL_RELATIONSHIP ) ImGui::Text( " Type: Parallel" );
        if ( corr->RelationshipType() == ANTIPARALLEL_RELATIONSHIP ) ImGui::Text( " Type: Antiparallel" );
        if ( corr->RelationshipType() == PERPENDICULAR_RELATIONSHIP ) ImGui::Text( " Type: Perpendicular" );

        FETFeature * f1 = corr->Feature(0);
        FETFeature * f2 = corr->Feature(1);
        Proxy * p1 = NULL; 
        Proxy * p2 = NULL;
        auto shape_proxy1 = seq_data.model->shape_to_proxy->find( f1->Shape() );
        auto shape_proxy2 = seq_data.model->shape_to_proxy->find( f2->Shape() );
        if( shape_proxy1 != seq_data.model->shape_to_proxy->end() )
        {
          p1 = shape_proxy1->second;
        }
        if( shape_proxy2 != seq_data.model->shape_to_proxy->end() )
        {
          p2 = shape_proxy2->second;
        }

        if( p1 && p2 )
        {
          ImGui::Text(" %d <-> %d", p1->structure_index, p2->structure_index );
          ImGui::Text(" Affinity: %f", corr->affinity );
          R3Vector n1 =  p1->feature->Normal();
          R3Vector n2 =  p2->feature->Normal();
          RNScalar angle = R3InteriorAngle( n1, n2 );
          ImGui::Text(" Angle: %f", angle );
        } 
        
      }
      else if( class_id == TRAJECTORY_ID )
      {
        ImGui::Text("Trajectory");
        double distance = seq_data.parametrization[view_opts.frame_idx_B] - 
                         seq_data.parametrization[view_opts.frame_idx_A];
        ImGui::Text(" Start: %d", view_opts.frame_idx_B );
        ImGui::Text(" End: %d", view_opts.frame_idx_A );
        ImGui::Text(" Extend: %f", distance );
      }
      else 
      {
        ImGui::Text("Error!");
      }
      ImGui::End();
    }
    else
    {
      if( lmb_clicked )
      {
        view_opts.selected_class_id = -1;
        view_opts.selected_object_id = -1;
        render.geometrical_corrs_dirty = true;
      }
    }
  }
  else
  {
    view_opts.selected_class_id = -1;
    view_opts.selected_object_id = -1;
  }

  
  // begin rendering interface
  ImGui::SetNextWindowPos(ImVec2(bsc::g_window.size.x - 256, 0));
  ImGui::SetNextWindowSize(ImVec2(256, bsc::g_window.size.y));

  ImGui::Begin("Options", NULL, ImVec2(-1, -1), -1.0f,
               ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                   ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);

  // Primitive Tabs
  const char *names[2] = {"Optimization", "View"};
  int tab_idx = gui::tabs(names, 2, 256);

  // Actual gui starts here
  if( tab_idx == 0)
  {
    if( ImGui::Button("Optimize", ImVec2(236, 24)))
    {
      run_optimization( &seq_data, &optim_params, &view_opts, &opts );
      rendering::request_update( &render, &opts );
    }

    if( ImGui::Button("Align with World", ImVec2(236, 24)))
    {
      align::reconstruction_to_world(&seq_data, &opts, 4);
      align::reconstruction_to_world(&seq_data, &opts, 4);
      update_bboxes( &seq_data, &opts, optim_params.max_depth );
      view_opts.xform_idx = seq_data.n_xforms - 1;
      view_opts.xforms_ptr    = &(seq_data.xforms_history[ 
                                  view_opts.xform_idx * seq_data.n_frames ] );
      rendering::request_update( &render, &opts );
    }

    if( ImGui::Button("Write Ply", ImVec2(236, 24)))
    {
      printf( "%s\n", opts.conf_filename );
      char buffer[ 128 ];
      char * dot_ptr = strrchr( opts.conf_filename, '.' );
      printf("%d\n", (int)(dot_ptr - opts.conf_filename));
      strncpy( buffer, opts.conf_filename, (int)(dot_ptr - opts.conf_filename));
      buffer[(int)(dot_ptr - opts.conf_filename)] = '\0';
      sprintf( buffer, "%s_%02d_%d_%d_%5.3f_%d_%d.ply", buffer,
                                       optim_params.iter_idx,
                                       (int)optim_params.cp_weight,
                                       (int)optim_params.traj_weight,
                                       optim_params.inertia_weight,
                                       (int)optim_params.geometrical_weight,
                                       (int)optim_params.hierarchical_weight );
      io::write_ply( buffer, seq_data.reconstruction, &opts, 4.5);
    }

    if( opts.output_conf_filename != NULL) 
    {
      if( ImGui::Button("Save Results", ImVec2(236, 24)))
      {
        io::write_configuration_file( opts.output_conf_filename, &seq_data, &opts );

        if( opts.xforms_history_filename )
        {
          io::write_transformation_history( opts.xforms_history_filename,
                                            &seq_data,
                                            &opts );
        }
        
        if( strcmp( seq_data.dataset, "tum" ) == 0 )
        {
          io::write_tum_trajectory_file("tum.traj", &seq_data, &opts );
        }

      }
    }

    if( ImGui::CollapsingHeader("Weights"))
    {
      ImGui::InputFloat( "Closest Point", &optim_params.cp_weight, 0.0, 10000.0 );
      ImGui::InputFloat( "Trajectory",    &optim_params.traj_weight, 0.0, 10000.0 );
      ImGui::InputFloat( "Inertia",       &optim_params.inertia_weight, 0.0, 10000.0 );
      ImGui::InputFloat( "Structure",     &optim_params.geometrical_weight, 0.0, 10000.0 );
      ImGui::InputFloat( "To Model",      &optim_params.hierarchical_weight, 0.0, 10000.0 );
    }

    // Optimization related options
    if( ImGui::CollapsingHeader("Info"))
    {
      gui::sequence_info(&seq_data);
    }


    if( ImGui::CollapsingHeader("Errors") )
    {
      if( seq_data.n_gt_corrs > 0 )
      {
        gui::optimization_errors( &seq_data );

        if( ImGui::CollapsingHeader("Show Correspondences") )
        {
          gui::correspondence_buttons( &view_opts, &seq_data );
        }
      }
    }
    if( ImGui::CollapsingHeader("Threshold function") )
    {
      static float max_dist = 0;
      static int shape_idx = 0;
      static int iter_idx = 0;
      float prev_max_dist = max_dist;
      int prev_shape_idx = shape_idx;
      int prev_iter_idx = iter_idx;
      ImGui::SliderFloat( "Max Dist", &max_dist, 0, 1000 );
      ImGui::SliderInt( "Shape Idx.", &shape_idx, 0, opts.frame_end );
      ImGui::SliderInt( "Iter Idx.", &iter_idx, 0, 10 );

      static std::vector<float> threshold_values( opts.frame_end, 0.5f );
      float min = std::min( optim_params.init_distance_threshold, optim_params.final_distance_threshold );
      float max = std::max( optim_params.init_distance_threshold, optim_params.final_distance_threshold );
      if( shape_idx != prev_shape_idx || max_dist != prev_max_dist  || iter_idx != prev_iter_idx ) {
        std::vector<double> parametrization( opts.frame_end );
        parametrization[0] = 0.0f;
        for( int i = 1; i < opts.frame_end; i++ )
        {
          bsc::vec3d origin0 = bsc::vec3d(seq_data.current_xforms[i-1][3]);
          bsc::vec3d origin1 = bsc::vec3d(seq_data.current_xforms[i][3]);
          RNLength d = bsc::norm(origin0 - origin1);
          parametrization[i] = parametrization[i - 1] + d;
        }

        int x = 0;
        for ( float &y : threshold_values )
        {
          y = optcorrs::calculate_threshold( shape_idx, x, iter_idx, max_dist, 
                                             parametrization, min, max );
          x++;
        }
      }
      ImGui::PlotLines("Threshold Function", threshold_values.data(), 
                                             threshold_values.size(),
                                             0, NULL, min, 2 * max, 
                                             ImVec2(256, 128));
    }
  }
  else
  {
    gui::view_options(&view_opts, &render, viewport, &seq_data, &opts);
  }
  ImGui::End();

  if( view_opts.show_fps )
  {
    gui::fps(viewport, &seq_data, &view_opts, total_n_points);
  }

  // Finish frame
  ImGui::Render();

  bsc::ui::end_frame();
}
#endif

////////////////////////////////////////////////////////////////////////////////
// Entry point
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{ 
  initialize_options( &opts, &view_opts , &optim_params );
  parse_arguments( argc, argv, &opts, &optim_params );
  initialize_random_numbers( &optim_params );


#ifdef FETREGISTER_USE_WINDOW
  if( !opts.batch_mode )
  {
    bsc::create_window("Camera Registration", 1280 + 256, 768 + 256,
                        bsc::vec3(0.47f, 0.49f, 0.53f),
                        bsc::vec3(0.2f, 0.21f, 0.25f) );
    bsc::set_init_funct(Init);
    bsc::set_display_funct(Display);
    return bsc::main_loop();
  }
  else 
#endif
  {
    if( initialize_optimization( &seq_data, &optim_params, &opts ) )
    {
      run_optimization( &seq_data, &optim_params, &view_opts, &opts );
      if( opts.output_conf_filename != NULL) 
      {
        io::write_configuration_file( opts.output_conf_filename, &seq_data, &opts );
      }
      if( opts.write_ply_model )
      {
        char buffer[ 128 ];
        char * dot_ptr = strrchr( opts.conf_filename, '.' );
        strncpy( buffer, opts.conf_filename, 
                                           (int)(dot_ptr - opts.conf_filename));
        buffer[(int)(dot_ptr - opts.conf_filename)] = '\0';
        sprintf( buffer, "%s_%02d_%d_%d_%5.3f_%d_%d.ply", buffer,
                                        optim_params.iter_idx,
                                        (int)optim_params.cp_weight,
                                        (int)optim_params.traj_weight,
                                        optim_params.inertia_weight,
                                        (int)optim_params.geometrical_weight,
                                        (int)optim_params.hierarchical_weight );
        io::write_ply( buffer, seq_data.reconstruction, &opts, 3.0);
      }
    }
    else
    {
      return 0;
    }
  }

  return 1;
}
