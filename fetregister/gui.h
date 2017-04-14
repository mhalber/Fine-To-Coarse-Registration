// FineToCoarseRegistration - fetregister - gui
//
// This file constains functions for generating user interface in windowed mode.
// It relies heavily on the 'dear, imgui' library by Omar Cornut.

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

#ifndef GUI_H_
#define GUI_H_

////////////////////////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////////////////////////

namespace gui
{
  i32 tabs( const char ** tab_names, const i32 n_tabs, const i32 width );
  i32 fps( const bsc::vec4i viewport, SequenceData * data, ViewOptions * view_opts, const i64 n_points = -1 );
  i32 sequence_info( const SequenceData * seq_data );
  i32 view_options( ViewOptions * view_opts, 
                    RenderingState * render,
                    bsc::vec4i viewport,
                    const SequenceData * seq_data,
                    const Options * opts );
  i32 correspondence_buttons( ViewOptions * view_opts, 
                              const SequenceData * seq_data );
  i32 optimization_errors( const SequenceData * seq_data );
}

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////
#ifdef GUI_IMPLEMENTATION

// primitive and very not-general tab support
i32 gui::
tabs( const char ** tab_names, const i32 n_tabs, const i32 width )
{
  ImDrawList *draw_list = ImGui::GetWindowDrawList();
  draw_list->PushClipRectFullScreen();
  static ImVec4 col = ImVec4(0.12f, 0.12f, 0.12f, 1.0f);
  const ImU32 col32 = ImColor(col);
  ImVec2 p( bsc::g_window.size.x - 256, 0 );
  draw_list->AddRectFilled(ImVec2(p.x, p.y), ImVec2(p.x + 256, p.y + 42), col32, 0.0f);
  draw_list->PopClipRect();

  ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 0.0f);

  static i32 active_tab_idx = 0; // this is the trick, this guy is persistent
  i32 prev_tab_idx = active_tab_idx;
  ImGuiStyle & stl = ImGui::GetStyle();
  ImVec2 tab_size;
  tab_size.x = (width - (stl.ItemSpacing.x * n_tabs - 1) - stl.WindowPadding.x * 2) / n_tabs;
  tab_size.y = 32;
  
  for ( int i = 0 ; i < n_tabs ; ++i )
  {
    if ( prev_tab_idx == i )
    {
      ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.2f, 0.2f, 0.2f, 1.00f));
      ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.2f, 0.2f, 0.2f, 1.00f));
      ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.2f, 0.2f, 0.2f, 1.00f));
    }
    if ( ImGui::Button( tab_names[i], tab_size ) )
    {
      active_tab_idx = i;
    }
    if ( prev_tab_idx == i )
    {
      ImGui::PopStyleColor(3);
    }
    if ( i != n_tabs - 1 ) ImGui::SameLine();
  }

  ImGui::PopStyleVar(1);
  ImGui::SetCursorPosY(ImGui::GetCursorPosY() + 5);

  return active_tab_idx;
}

i32 gui::
sequence_info( const SequenceData * seq_data )
{
  ImGui::Text("Format           : %s", seq_data->format );
  ImGui::Text("Dataset          : %s", seq_data->dataset );
  ImGui::Text("No. of Frames    : %d\n", seq_data->n_frames );
  ImGui::Text("Color Res        : %d %d", seq_data->color_res.x, 
                                          seq_data->color_res.y);
  ImGui::Text("Depth Res        : %d %d", seq_data->depth_res.x, 
                                          seq_data->depth_res.y);

  bsc::mat3 cK = seq_data->color_intrinsics;
  ImGui::Text("Color Intrinsics : "); 
  ImGui::Text("\t%8.2f %8.2f %8.2f", cK[0][0], cK[0][1], cK[0][2]);
  ImGui::Text("\t%8.2f %8.2f %8.2f", cK[1][0], cK[1][1], cK[1][2]);
  ImGui::Text("\t%8.2f %8.2f %8.2f", cK[2][0], cK[2][1], cK[2][2]);

  bsc::mat3 dK = seq_data->depth_intrinsics;
  ImGui::Text("Depth Intrinsics :"); 
  ImGui::Text("\t%8.2f %8.2f %8.2f", dK[0][0], dK[0][1], dK[0][2]);
  ImGui::Text("\t%8.2f %8.2f %8.2f", dK[1][0], dK[1][1], dK[1][2]);
  ImGui::Text("\t%8.2f %8.2f %8.2f", dK[2][0], dK[2][1], dK[2][2]);

  bsc::mat4 cT = seq_data->color_extrinsics;
  ImGui::Text("Color Extrinsics :"); 
  ImGui::Text("\t%8.2f %8.2f %8.2f %8.2f", cT[0][0], cT[0][1], cT[0][2], cT[0][3]);
  ImGui::Text("\t%8.2f %8.2f %8.2f %8.2f", cT[1][0], cT[1][1], cT[1][2], cT[1][3]);
  ImGui::Text("\t%8.2f %8.2f %8.2f %8.2f", cT[2][0], cT[2][1], cT[2][2], cT[2][3]);
  ImGui::Text("\t%8.2f %8.2f %8.2f %8.2f", cT[3][0], cT[3][1], cT[3][2], cT[3][3]);

  bsc::mat4 dT = seq_data->depth_extrinsics;
  ImGui::Text("Depth Extrinsics :"); 
  ImGui::Text("\t%8.2f %8.2f %8.2f %8.2f", dT[0][0], dT[0][1], dT[0][2], dT[0][3]);
  ImGui::Text("\t%8.2f %8.2f %8.2f %8.2f", dT[1][0], dT[1][1], dT[1][2], dT[1][3]);
  ImGui::Text("\t%8.2f %8.2f %8.2f %8.2f", dT[2][0], dT[2][1], dT[2][2], dT[2][3]);
  ImGui::Text("\t%8.2f %8.2f %8.2f %8.2f", dT[3][0], dT[3][1], dT[3][2], dT[3][3]);

  return 1;
}

i32 gui::
view_options( ViewOptions * view_opts,
              RenderingState * render,
              bsc::vec4i viewport,
              const SequenceData * seq_data,
              const Options * opts )
{
  if (ImGui::CollapsingHeader("View Options",
                              ImGuiTreeNodeFlags_DefaultOpen))
  {

    ImGui::Checkbox("Show Frames",     &view_opts->show_frames);
    if( view_opts->show_frames )
    {
      ImGui::InputInt("Frame Skip", &view_opts->frame_skip, 5, 5 );
    }
    ImGui::Checkbox("Show Cameras",    &view_opts->show_cameras);
    ImGui::Checkbox("Show Trajectory", &view_opts->show_trajectory);
    ImGui::Checkbox("Show Bboxes",     &view_opts->show_bboxes);
    ImGui::Checkbox("Inspection Mode", &view_opts->inspection_mode);
    ImGui::SliderInt("Feat. Disp. Mode", (int*)&view_opts->feat_disp_mode, 0, 4 );
    ImGui::SliderFloat("Point Size", &view_opts->point_size, 1.0f, 100.0f);
    ImGui::SliderFloat("Max Dist",   &view_opts->max_dist, 0.1f, 10.0f);
    ImGui::SliderFloat("Max Height", &view_opts->max_height, -15.0f, 15.0f);

    ImGui::Separator();

    ImGui::Checkbox("Show Axes", &view_opts->show_axes);
    ImGui::Checkbox("Show Grid", &view_opts->show_grid);
    ImGui::Checkbox("Show FPS",  &view_opts->show_fps);

    ImGui::Separator();

    r32 prev_offset = view_opts->proxy_offset;
    ImGui::Checkbox("Show C.P. Corrs.",           &view_opts->show_cp_corrs);
    ImGui::Checkbox("Show Hierarchical. Corrs.",  &view_opts->show_hierarchical_corrs);
    ImGui::Checkbox("Show Geometrical. Corrs.",   &view_opts->show_geometrical_corrs);
    ImGui::Checkbox("Show G.T. Corrs.",           &view_opts->show_gt_corrs);
    ImGui::Checkbox("Show Proxies",               &view_opts->show_proxies);
    if ( view_opts->show_proxies )
    {
      ImGui::SliderInt("Proxy Level", &view_opts->proxy_level, -1, 1 );
    }
    if ( view_opts->show_hierarchical_corrs )
    {
      if( ImGui::CollapsingHeader("Show Correspondences Type") )
      {
        ImGui::Checkbox("Coincident", &view_opts->show_relationships_by_type[COINCIDENT_RELATIONSHIP]);
        ImGui::Checkbox("Parallel", &view_opts->show_relationships_by_type[PARALLEL_RELATIONSHIP]);
        ImGui::Checkbox("Antiparallel", &view_opts->show_relationships_by_type[ANTIPARALLEL_RELATIONSHIP]);
        ImGui::Checkbox("Perpendicular", &view_opts->show_relationships_by_type[PERPENDICULAR_RELATIONSHIP]);
      }
    }
    ImGui::SliderInt("Proxy. Disp. Mode", (int*)&view_opts->proxy_disp_mode, 0, 2 );
    ImGui::SliderFloat("Offset", &view_opts->proxy_offset, 0.0f, 2.0f);
    if ( ( view_opts->show_geometrical_corrs || view_opts->show_hierarchical_corrs ) &&
         prev_offset != view_opts->proxy_offset )
    {
      render->geometrical_corrs_dirty = true;
      render->hierarchical_corrs_dirty= true;
    }
    ImGui::SliderFloat("Line Width", &view_opts->line_width, 1.0f, 20.0f);
    
    ImGui::Separator();
    ImGui::Text("Frame Index Limits");
    i32 prev_idx_A = view_opts->frame_idx_A;
    i32 prev_idx_B = view_opts->frame_idx_B;
    i32 prev_show_only_pair = view_opts->show_only_pair;
    ImGui::SliderInt( "Min. Idx", &view_opts->frame_idx_A,
                        0, opts->frame_end-2 );
    ImGui::PushButtonRepeat( true );
    if ( ImGui::Button("<<##idx_A", ImVec2( 120, 16 ) ) ) { 
      view_opts->frame_idx_A--;
    }
    ImGui::SameLine();
    if ( ImGui::Button(">>##idx_A", ImVec2( 120, 16 ) ) ) { 
      view_opts->frame_idx_A++;
    }
    ImGui::PopButtonRepeat();    
    ImGui::SliderInt( "Max. Idx", &view_opts->frame_idx_B,
                        0, opts->frame_end-1 );

    ImGui::PushButtonRepeat( true );  
    if ( ImGui::Button("<<##idx_B", ImVec2( 120, 16 ) ) ) { 
      view_opts->frame_idx_B--; 
    } 
    ImGui::SameLine();
    if ( ImGui::Button(">>##idx_B", ImVec2( 120, 16 ) ) ) { 
      view_opts->frame_idx_B++; 
    }
    ImGui::PopButtonRepeat();    
    ImGui::Checkbox( "Show Only Pair ", &view_opts->show_only_pair );
    if ( view_opts->frame_skip > 1 )
    {
      view_opts->frame_idx_A = view_opts->frame_idx_A -
                              (view_opts->frame_idx_A % view_opts->frame_skip);
      view_opts->frame_idx_B = view_opts->frame_idx_B -
                              (view_opts->frame_idx_B % view_opts->frame_skip);
    }
    view_opts->frame_idx_A = std::max( std::min(view_opts->frame_idx_A,
                                                     opts->frame_end-2 ), 0 );
    view_opts->frame_idx_B = std::max( std::min(view_opts->frame_idx_B,
                                                     opts->frame_end-1 ), 0 );
    if ( prev_idx_A != view_opts->frame_idx_A ||
         prev_idx_B != view_opts->frame_idx_B ||
         prev_show_only_pair != view_opts->show_only_pair )
    {
      render->gt_corrs_dirty             = true;
      // render->cp_corrs_dirty             = true;
      // render->loop_closure_matches_dirty = true;
      // render->geometrical_corrs_dirty       = true;
      // render->hierarchical_corrs_dirty        = true;
      render->trajectory_dirty           = true;
    }
    ImGui::Separator();
  
    i32 prev_idx = view_opts->xform_idx;
    ImGui::SliderInt("Xform Idx", &view_opts->xform_idx, 
                                                    0, seq_data->n_xforms - 1 );
    if ( ImGui::Button("<<##idx_xform", ImVec2( 120, 16 ) ) ) { 
      view_opts->xform_idx = std::max( 0, view_opts->xform_idx - 1); 
    }
    ImGui::SameLine();
    if ( ImGui::Button(">>##idx_xform", ImVec2( 120, 16 ) ) ) { 
      view_opts->xform_idx = std::min ( seq_data->n_xforms-1, 
                                                     view_opts->xform_idx + 1 );
    }
    if ( prev_idx != view_opts->xform_idx )
    {
      view_opts->xforms_ptr = &(seq_data->xforms_history[ 
                                  view_opts->xform_idx * seq_data->n_frames ] );
      render->gt_corrs_dirty             = true;
      // render->cp_corrs_dirty             = true;
      // render->loop_closure_matches_dirty = true;
      // render->geometrical_corrs_dirty       = true;
      // render->hierarchical_corrs_dirty        = true;
      render->trajectory_dirty           = true;
    }
  
    if ( ImGui::CollapsingHeader("Other Methods") )
    {
      bool clicked = false;
      if ( seq_data->rgbdsfm_xforms )
      {
        if ( ImGui::Button("Sun3d", ImVec2(105, 32)) )
        {
          view_opts->xforms_ptr = seq_data->rgbdsfm_xforms;
          clicked = true;
        }
        ImGui::SameLine();
      }
      if ( seq_data->robust_reconstruction_xforms )
      {
        if ( ImGui::Button("Robust\n Reconstr.", ImVec2(105, 32)) )
        {
          view_opts->xforms_ptr = seq_data->robust_reconstruction_xforms;
          clicked = true;
        }
      }

      if ( seq_data->kintinuous_xforms )
      {
        if ( ImGui::Button("Kintinuous", ImVec2(105, 32)) )
        {
          view_opts->xforms_ptr = seq_data->kintinuous_xforms;
          clicked = true;
        }
        ImGui::SameLine();
      }

      if ( seq_data->elastic_fusion_xforms )
      {
        if ( ImGui::Button("Elastic\n Fusion", ImVec2(105, 32)) )
        {
          view_opts->xforms_ptr = seq_data->elastic_fusion_xforms;
          clicked = true;
        }
      }
      
      if ( ImGui::Button("Ours", ImVec2(215, 32)) )
      {
        view_opts->xforms_ptr = seq_data->current_xforms; 
        view_opts->xform_idx = seq_data->n_xforms - 1;
        
        clicked = true;
      }
      if ( seq_data->gt_xforms )
      {
        if ( ImGui::Button("Ground Truth", ImVec2(215, 32)) )
        {
          view_opts->xforms_ptr = seq_data->gt_xforms; 
          clicked = true;
        }
      }

      if ( clicked )
      {
        render->trajectory_dirty     = true;
        render->gt_corrs_dirty       = true;
        render->cp_corrs_dirty       = true;
        render->geometrical_corrs_dirty = true;
        render->hierarchical_corrs_dirty  = true;
      }
    }
  }

  if (ImGui::CollapsingHeader("Camera"))
  {
    ImGui::SliderFloat("Vertical FOV", &view_opts->fov_y, 0.1f, 120.0f);
    if (ImGui::Button("Front", ImVec2(74, 24)))
    {
      render->cam_controls.set_front_view(&render->cam, &render->view);
    }
    ImGui::SameLine();
    if (ImGui::Button("Right", ImVec2(74, 24)))
    {
      render->cam_controls.set_right_view(&render->cam, &render->view);
    }
    ImGui::SameLine();
    if (ImGui::Button("Top", ImVec2(74, 24)))
    {
      render->cam_controls.set_top_view(&render->cam, &render->view);
    }

    if (ImGui::Button("Back", ImVec2(74, 24)))
    {
      render->cam_controls.set_back_view(&render->cam, &render->view);
    }
    ImGui::SameLine();
    if (ImGui::Button("Left", ImVec2(74, 24)))
    {
      render->cam_controls.set_left_view(&render->cam, &render->view);
    }
    ImGui::SameLine();
    if (ImGui::Button("Bottom", ImVec2(74, 24)))
    {
      render->cam_controls.set_bottom_view(&render->cam, &render->view);
    }
  }
  if (ImGui::CollapsingHeader("Background"))
  {
    ImGui::ColorEdit3("Gradient Top", &(bsc::g_window.bckgrd_col_top.r));
    ImGui::ColorEdit3("Gradient Bottom", &(bsc::g_window.bckgrd_col_bot.r));
  }
  if ( ImGui::Button( "Take Screenshot", ImVec2( 236, 24 ) ) )
  {
    bsc::take_screenshot( viewport );
  }
  if ( ImGui::Button( "Save Ply", ImVec2( 236, 24 ) ) )
  {
    io::write_ply( "test.ply", seq_data, opts, view_opts->max_dist );
  }
  return 1;
}

i32 gui::
fps( const bsc::vec4i viewport, SequenceData * data, ViewOptions * view_opts, const i64 n_points )
{
  char n_points_string[512];
  char image_name_A[512];
  char image_name_B[512];
  sprintf( n_points_string, "NPoints  : %lld", n_points );
  sprintf( image_name_A, "%s", data->depth_names[view_opts->frame_idx_A] );
  sprintf( image_name_B, "%s", data->depth_names[view_opts->frame_idx_B] );
  ImVec2 add_size( std::max( ImGui::CalcTextSize( n_points_string ).x,
                             ImGui::CalcTextSize( image_name_A ).x),
                   ImGui::CalcTextSize( n_points_string ).y +
                   ImGui::CalcTextSize( image_name_A ).y +
                   ImGui::CalcTextSize( image_name_B ).y);
  ImVec2 size( 132, 28 );
  if ( n_points > 0 )
  {
    size.x = std::max( add_size.x + 16.0f, size.x );
    size.y += add_size.y + 8.0f;
  }

  ImGui::SetNextWindowPos( ImVec2( viewport.z - size.x,
                                   viewport.w - size.y ) );
  
  ImGui::Begin("StatusWindow", NULL, ImVec2(-1, -1), 0.0f,
                ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                    ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings);
  ImGui::Text("%5.4fms (%dFPS)\n", 1.0 / ImGui::GetIO().Framerate,
              (i32)ImGui::GetIO().Framerate);
  ImGui::Text("%s", image_name_A);
  ImGui::Text("%s", image_name_B);
  
  if ( n_points > 0 )
  {
    ImGui::Text("%s", n_points_string );
  }

  ImGui::End();
  return 1;
}

i32 gui::
correspondence_buttons( ViewOptions * view_opts, const SequenceData * seq_data )
{
  // SORT CORRESPONDENCES SOMEHOW -> ONE POSSIBILITY ARE TUPLES

  // ok this is a  workable solution --> I will probably need to do some benchmarks
  // Sort function should return sorted array and index array

  tuple * tuples = (tuple*)malloc( sizeof(tuple) * seq_data->n_gt_corrs );
  for ( int i = 0 ; i < seq_data->n_gt_corrs ; ++i )
  {
    tuple t;
    t.val = seq_data->per_correspondence_error[i];
    t.idx = i;
    tuples[i] = t;
  }
  qsort( &(tuples[0]), seq_data->n_gt_corrs, sizeof(tuple), idx_compare );

  for ( i32 idx = 0; 
        idx < seq_data->n_gt_corrs ; 
        ++idx )
  {
    int corr_idx = tuples[ idx ].idx;

    const Correspondence * corr = &(seq_data->gt_corrs[corr_idx]);
    char msg[MAX_STR_LENGTH];
    sprintf( msg, "%d:%d-%d", corr_idx, corr->idx_A, corr->idx_B );
    if ( ImGui::Button(msg, ImVec2(70,24) ) )
    {
      view_opts->frame_idx_A = corr->idx_A;
      view_opts->frame_idx_B = corr->idx_B;
    }
    if ( idx == 0 || idx % 3 != 0 ) ImGui::SameLine();
  }

  free( tuples );
  return 1;
}

i32 gui::
optimization_errors( const SequenceData * seq_data )
{

  ImGui::Text("Current Error: %f", seq_data->errors[seq_data->n_xforms - 1]);
  
  r32 max_error = 0.0;
  for (i32 i = 0; i < seq_data->n_xforms; ++i)
  {
    if (seq_data->errors[i] > max_error)
      max_error = seq_data->errors[i];
  }
  ImGui::PlotLines( "Errors", &(seq_data->errors[0]), seq_data->n_xforms, 
                    0, NULL, 0.0, max_error, ImVec2(236, 128) );
  ImGui::PlotHistogram("PerCorrespondence", 
                      &(seq_data->per_correspondence_error[0]), 
                        seq_data->n_gt_corrs,
                          0, NULL, 0.0, 0.2, ImVec2( 236, 128 ) );
  ImGui::Separator();

  ImGui::Columns(4, "errors"); 
  ImGui::Separator();
  ImGui::Text("Method"); ImGui::NextColumn();
  ImGui::Text("RMSE"); ImGui::NextColumn();
  ImGui::Text("Mean"); ImGui::NextColumn();
  ImGui::Text("Std. Dev"); ImGui::NextColumn();
  ImGui::Separator();
  // calculate errors
  const char* method_names[5] = { "Ours", "Sun3D", "Robust Reconstruction", 
                                  "Elastic Fusion", "Kintinuous" };
  r32 rmse[5], mean[5], sdev[5];
  eval::calculate_error_stats( seq_data, &(rmse[0]), &(mean[0]), &(sdev[0]) );

  // find minimum
  i32 min_rmse_idx = -1, min_mean_idx = -1, min_sdev_idx = -1;
  r32 min_rmse = 1e9, min_mean = 1e9, min_sdev = 1e9;
  for ( i32 i = 0 ; i < 5; ++i )
  {
    if ( rmse[i] < min_rmse ) { min_rmse = rmse[i], min_rmse_idx = i; }
    if ( mean[i] < min_mean ) { min_mean = mean[i], min_mean_idx = i; }
    if ( sdev[i] < min_sdev ) { min_sdev = sdev[i], min_sdev_idx = i; }
  }

  // draw table
  ImVec4 default_color = ImVec4( 0.9f, 0.9f, 0.9f, 1.0f );
  ImVec4 best_color    = ImVec4( 1.0f, 0.5f, 0.0f, 1.0f );
  for ( i32 i = 0; i < 5; ++i )
  {
    ImGui::Text("%s\n", method_names[i]); ImGui::NextColumn();
    
    ImVec4 rmse_color = default_color;
    if ( i == min_rmse_idx ) rmse_color = best_color;
    ImGui::TextColored(rmse_color, "%1.4f", rmse[i]); ImGui::NextColumn();
    
    ImVec4 mean_color = default_color;
    if ( i == min_mean_idx ) mean_color = best_color;
    ImGui::TextColored(mean_color, "%1.4f", mean[i]); ImGui::NextColumn();
   
    ImVec4 sdev_color = default_color;
    if ( i == min_sdev_idx ) sdev_color = best_color;
    ImGui::TextColored(sdev_color, "%1.4f", sdev[i]); ImGui::NextColumn();
  }
  ImGui::Columns(1);
  ImGui::Separator();
  return 1;
}


#endif // GUI_IMPlEMENTATION

#endif //GUI_H_