// FineToCoarseRegistration - fetregister - rendering
//
// This file contains routines for drawing the pointcloud and interface when
// running fetregister in windowed mode.



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


#ifndef RENDERING_H_
#define RENDERING_H_

////////////////////////////////////////////////////////////////////////////////
// Interface / Global State
////////////////////////////////////////////////////////////////////////////////

enum ClassIDs
{
  POINTCLOUD_ID = 1,
  PROXY_ID,
  CP_CORR_ID,
  GEOMETRICAL_CORR_ID,
  HIERARCHICAL_CORR_ID,
  TRAJECTORY_ID,
};

// TODO: Clean up arguments
struct RenderingState
{
  // view control
  bsc::camera cam;
  bsc::trackball_controls cam_controls;
  bsc::mat4 view;
  bsc::mat4 proj;

  // geometries we will draw
  bsc::gpu_geometry * pointclouds_geo = NULL;
  bsc::gpu_geometry * bboxes_geo      = NULL;
  bsc::gpu_geometry gt_corrs_geo;
  bsc::gpu_geometry cp_corrs_geo;
  bsc::gpu_geometry geometrical_corrs_geo;
  bsc::gpu_geometry hierarchical_corrs_geo;
  bsc::gpu_geometry trajectory_geo;
  bsc::gpu_geometry axes_geo;
  bsc::gpu_geometry frustum_geo;
  bsc::gpu_geometry grid_geo;
  bsc::gpu_geometry disc_geo;

  // shaders
  bsc::shader_prog default_shader;
  bsc::shader_prog line_adj_shader;
  bsc::shader_prog line_shader;
  bsc::shader_prog pointcloud_shader;
  bsc::shader_prog fullscreen_shader;

  // state flags
  bool * pointclouds_dirty        = NULL;
  bool * bboxes_dirty             = NULL;
  bool gt_corrs_dirty             = true;
  bool cp_corrs_dirty             = true;
  bool hierarchical_corrs_dirty         = true;
  bool geometrical_corrs_dirty       = true;
  bool trajectory_dirty           = true;

  // framebuffers / renderbuffers
  bsc::vec4i viewport;
  float pixel_ratio;
  GLuint color_buffer;
  GLuint index_buffer;
  GLuint depth_renderbuffer;
  GLuint fullscreen_framebuffer;
  GLuint fullscreen_vao;
};

namespace rendering
{
  void init( RenderingState * render,
             const SequenceData * seq_data,
             const Options * opts,
             const bsc::vec4i viewport );
  void draw( RenderingState * render, 
             const SequenceData * seq_data,
             const ViewOptions * view_opts,
             const OptimizationParameters * optim_params,
             const Options * opts,
             const bsc::vec4i viewport );
  void request_update( RenderingState * render,
                       const Options * opts );

// TODO: Decide if init is a correct name even.
  void init_axes( bsc::gpu_geometry *cmd );

  void init_frustum( bsc::gpu_geometry *cmd, 
                     const bsc::mat3 *K );

  void init_grid( bsc::gpu_geometry *cmd, 
                  i32 grid_size, 
                  r32 grid_spacing = 1.0 );

  void init_disc( bsc::gpu_geometry *gpu_geo, 
                  i32 divs = 32, 
                  r32 redius = 0.25f );

  void init_trajectory( bsc::gpu_geometry *traj_geo,
                        const bsc::mat4d *extrinsics,
                        const int *n_inliers,
                        const Options *opts,
                        const i32 low_lim, 
                        const i32 high_lim );

  void init_framebuffer( RenderingState * render, bsc::vec4i viewport );

  i32 init_pointclouds_from_frames( RenderingState * render,
                                    const SequenceData * seq_data,
                                    const Options * opts );
  i32 init_pointclouds_from_reconstruction( RenderingState * render,
                                            const SequenceData * seq_data,
                                            const ViewOptions * view_opts,
                                            const Options * opts );


  void init_gt_corrs( bsc::gpu_geometry * geo,
                      const Correspondence * correspondences,
                      const i32 n_correspondences,
                      const bsc::mat4d * xforms,
                      const i32 lim_low = 0, 
                      const i32 lim_high = 1e9 );

  // FET DEPENDENCE
  void init_bboxes( bsc::gpu_geometry *geo,
                    bool *dirty,
                    const SequenceData *seq_data,
                    const Options *opts );

  void init_cp_corrs( bsc::gpu_geometry *geo,
                      const SequenceData *seq_data,
                      const ViewOptions *view_opts );

  void init_geometrical_corrs( bsc::gpu_geometry *geo,
                            const SequenceData *seq_data,
                            const ViewOptions *view_opts );

  void init_hierarchical_corrs( bsc::gpu_geometry *geo,
                          const SequenceData *seq_data,
                          const ViewOptions *view_opts );

  // helpers for processing frames
  void frame_to_pointcloud( bsc::geometry_data * pointcloud,
                            const bsc::img *color_img,
                            const bsc::img_u16 *depth_img,
                            const bsc::mat3 K,
                            const i32 downsample_factor = 2,
                            const bool bitshift = true,
                            const r32 near_dist = 0.01f,
                            const r32 far_dist = 10.0f );
  i32 init_pointcloud( bsc::gpu_geometry *cmd,
                       const bsc::img *color_img,
                       const bsc::img_u16 *depth_img,
                       const bsc::mat3 K,
                       const r32 downsample_factor = 2,
                       const bool bitshift = true,
                       const r32 near_dist = 0.01f,
                       const r32 far_dist = 10.0f  );
}

#ifdef RENDERING_IMPLEMENTATION
////////////////////////////////////////////////////////////////////////////////
// Colors
////////////////////////////////////////////////////////////////////////////////

// color maps from ColorBrewer.org
bsc::vec3i divergent_map[15] = {bsc::vec3i(213, 62, 79),
                               bsc::vec3i(204, 23, 79),
                               bsc::vec3i(244, 109, 67),
                               bsc::vec3i(253, 174, 97),
                               bsc::vec3i(254, 224, 139),
                               bsc::vec3i(255, 255, 191),
                               bsc::vec3i(230, 245, 152),
                               bsc::vec3i(171, 221, 164),
                               bsc::vec3i(102, 194, 165),
                               bsc::vec3i(50, 136, 189),
                               bsc::vec3i(25, 123, 219),
                               bsc::vec3i(15, 163, 239),
                               bsc::vec3i(68, 83, 129),
                               bsc::vec3i(112, 183, 179),
                               bsc::vec3i(221, 183, 29) };

bsc::vec3i qual_map[9] = {bsc::vec3i(166, 206, 227),
                          bsc::vec3i(31, 120, 180),
                          bsc::vec3i(178, 223, 138),
                          bsc::vec3i(51, 160, 44),
                          bsc::vec3i(251, 154, 153),
                          bsc::vec3i(227, 26, 28),
                          bsc::vec3i(253, 191, 111),
                          bsc::vec3i(255, 127, 0),
                          bsc::vec3i(202, 178, 214)};

bsc::vec3i transition_map[9] = { bsc::vec3i(255, 128, 0),
                                 bsc::vec3i(240, 107, 51),
                                 bsc::vec3i(226, 86, 102),
                                 bsc::vec3i(195, 71, 134),
                                 bsc::vec3i(165, 56, 165),
                                 bsc::vec3i(130, 70, 197),
                                 bsc::vec3i(96, 83, 228),
                                 bsc::vec3i(48, 105, 242),
                                 bsc::vec3i(0, 128, 255)};

bsc::vec3i rel_type_map[9] = {bsc::vec3i(255, 58, 17),
                              bsc::vec3i(178, 81, 62),
                              bsc::vec3i(34, 90, 204),
                              bsc::vec3i(184, 255, 43) };

bsc::vec3i red2blue[3] = { bsc::vec3i( 220, 44, 44 ),
                           bsc::vec3i( 142, 33, 218 ),
                           bsc::vec3i( 44, 40, 231) };

////////////////////////////////////////////////////////////////////////////////
// Helpers
////////////////////////////////////////////////////////////////////////////////

bool 
is_within( const i32 idx, const ViewOptions * view_opts )
{
  return (idx >= view_opts->frame_idx_A && idx <= view_opts->frame_idx_B);
}

bool
is_within( const Proxy * proxy, const ViewOptions * view_opts )
{
  std::vector<i32> shapes_ids;
  optproxy::get_shape_ind_connected_to_proxy( proxy, shapes_ids );
  i32 should_show = 0; 
  for ( i32 idx = 0;
        idx < shapes_ids.size();
        idx++ )
  {
    if ( is_within( shapes_ids[idx], view_opts ) )
    {
      should_show = 1;
      break;
    }
  }
  return should_show;
}

////////////////////////////////////////////////////////////////////////////////
// Shaders
////////////////////////////////////////////////////////////////////////////////

static const char* vs_fullscreen = (char *)SHADER_HEAD STR
(
  out vec2 v_uv;
  void main()
  {
    uint idx = uint( gl_VertexID  );
    gl_Position = vec4(
        (float( idx &  1U ) ) * 4.0 - 1.0,
        (float( idx >> 1U  ) ) * 4.0 - 1.0,
        0.0, 1.0);
    v_uv = gl_Position.xy * 0.5 + 0.5;
  }
);

static const char* fs_fullscreen = (char *)SHADER_HEAD STR
(
  uniform sampler2D u_color;
  // uniform sampler2D u_index;
  
  in vec2 v_uv;
  out vec4 frag_color;
  void main()
  {
    frag_color = texture( u_color, v_uv );
  }
);

char *vs_pointcloud = (char *)SHADER_HEAD STR(
    layout(location = 0) in vec3 position;
    layout(location = 1) in vec3 normal;
    layout(location = 4) in vec4 color_a;
    layout(location = 5) in vec4 color_b;
    uniform mat4 mv;
    uniform mat4 m;
    uniform mat4 p;
    uniform float max_height;
    uniform float point_size;
    uniform float max_dist;
    uniform int coloring_mode ;
    out vec4 v_color;

    void main() {
      vec4 view_pos = mv * vec4(position, 1.0);

      // if greater than max dist, simply reject
      if ( position.z < -max_dist ||
           (m*vec4(position, 1.0)).y > max_height )
      { 
        view_pos = vec4(1000.0, 1000.0, 1000.0, 1.0);
      }

      gl_Position = p * view_pos;
      gl_PointSize = point_size / (-view_pos.z);
      
      // coloring_mode == 0 // DEFAULT, just show colors
      v_color = vec4(color_a.rgb, 1.0f);

      if ( coloring_mode == 1 ) // GENERATOR TYPE
      {
        switch ( int( color_a.a * 255 ) )
        {
          case 1:  // SIFT
            v_color =  vec4(31/255.0f, 120/255.0f, 180/255.0f, 1.0f ); break;
          case 12: // SILHOUETTE
            v_color = vec4(51/255.0f, 160/255.0f, 44/255.0f, 1.0f ); break;
          case 15: // RIDGE
            v_color = vec4(253/255.0f, 191/255.0f, 111/255.0f, 1.0f ); break;
          case 16: // VALLEY
            v_color = vec4(255/255.0f, 127/255.0f, 0/255.0f, 1.0f ); break;
          case 22: // PLANE
            v_color = vec4(251/255.0f, 154/255.0f, 153/255.0f, 1.0f ); break;
          case 31: // STRUCTURE
            v_color = vec4(121/255.0f, 237/255.0f, 26/255.0f, 1.0f ); break;
          case 11 : // BORDER
            v_color = vec4(1.0, 0.0f, 0.0f, 1.0f ); break;
          case 13 : // SHADOW
            v_color = vec4( 0.0f, 0.0f, 1.0f, 1.0f ); break;
          default:
            v_color = vec4( 0.0f, 0.0f, 0.0f, 1.0f ); break;
        }
      }
      else if ( coloring_mode == 2 ) // PRIMITVE_IDX
      {
        v_color = vec4(color_b.rgb, 1.0f);
      }
      else if( coloring_mode == 3 ) // NORMALS
      {
        vec4 n = m * vec4( normal, 0.0f );
        n = normalize(n);
        v_color = vec4( 0.5f * (n.xyz + 1.0f) , 1.0f );
      }
      else if ( coloring_mode == 4 ) // PHONG
      {
        vec4 n = mv * vec4(normal, 0.0f );
        n = normalize(n);
        vec3 l = -view_pos.xyz;
        l = normalize(l);
        float ndotl = dot( vec3(n), l );
        v_color = vec4(ndotl, ndotl, ndotl, 1.0f);
        // cull backfacing points
        if (ndotl <= 0.0) 
        {
          view_pos = vec4(1000.0f, 1000.0f, 1000.0f, 1.0f);
          gl_Position = p * view_pos;
        }
      }
    });

char *fs_pointcloud = (char *)SHADER_HEAD STR(
    in vec4 v_color;
    layout (location = 0) out vec4 frag_color;
    layout (location = 1) out vec4 frag_id;
    uniform int class_id;
    uniform int object_id;
    void main() 
    {
      // cute hack of round points
        vec2 coord = gl_PointCoord - vec2(0.5);  //from [0,1] to [-0.5,0.5]
        if(length(coord) > 0.5)                  //outside of circle radius?
          discard;
        frag_color = v_color;
        frag_id = vec4( class_id, object_id, 0.0, 1.0);
    });

char *vs_default = (char *)SHADER_HEAD STR(
    layout(location = 0) in vec3 position;
    layout(location = 4) in vec4 color;
    uniform mat4 mvp;
    out vec4 v_color;
    void main() {
      gl_Position = mvp * vec4(position, 1.0);
      v_color = color;
    });

char *fs_default = (char *)SHADER_HEAD STR(
    in vec4 v_color;
    layout (location = 0) out vec4 frag_color;
    layout (location = 1) out vec4 frag_id;

    uniform bool use_uniform_color;
    uniform vec3 uniform_color;
    uniform int class_id;
    uniform int object_id;

    void main() {
      if (!use_uniform_color)
      {
        frag_color = v_color;
      }
      else
      {
        frag_color = vec4(uniform_color, 1.0f);
      }
      frag_id = vec4(class_id, object_id, 0.0, 1.0);
    });


// Shader program for drawing lines, as modern opengl dropped line width
// Adapted from
//https://forum.libcinder.org/topic/smooth-thick-lines-using-geometry-shader
// Fixed issue with clipping

char * vs_lines = (char *) SHADER_HEAD STR(
    layout(location = 0) in vec3 position;
    layout(location = 4) in vec4 color;
    layout(location = 8) in float user_data;
    uniform mat4 mvp;
    out vec4 v_color;
    flat out float v_object_id;
    void main() 
    {
      gl_Position = mvp * vec4(position, 1.0);
      v_color = color;
      v_object_id = user_data;
    });

char * gs_line_adj = (char *) SHADER_HEAD STR(
  layout( lines_adjacency ) in;
  layout( triangle_strip, max_vertices = 10 ) out;

  in vec4 v_color[];
  flat in float v_object_id[];
  out vec4 g_color;
  flat out float g_object_id;

  uniform float	line_width;		// the thickness of the line in pixels
  uniform float	miter_limit;	// 1.0: always miter, -1.0: never miter, 0.75: default
  uniform vec2	win_size;		  // the size of the viewport in pixels

  vec2 to_screen_space(vec4 pos)
  {
    // Division by w takes us to normalized device coordinates and then
    // win scale is simply the viewport transform:
    // http://www.songho.ca/opengl/gl_transform.html
    return vec2( pos.xy / pos.w ) * win_size;
  }

  vec4 to_clip_space( vec4 pos )
  {
    return vec4( (pos.xy / win_size) * pos.w , pos.zw );
  }

  void main(void)
  {
    // get the four vertices passed to the shader:
    vec4 p0 = gl_in[0].gl_Position; 
    vec4 p1 = gl_in[1].gl_Position; 
    vec4 p2 = gl_in[2].gl_Position; 
    vec4 p3 = gl_in[3].gl_Position; 

    //move to screen space
    vec2 sp0 = to_screen_space( p0 );	// start of previous segment
    vec2 sp1 = to_screen_space( p1 );	// end of previous segment, start of current segment
    vec2 sp2 = to_screen_space( p2 );	// end of current segment, start of next segment
    vec2 sp3 = to_screen_space( p3 );	// end of next segment
    
    // determine the direction of each of the 3 segments (previous, current, next)
    vec2 v0 = normalize( sp1 - sp0 );
    vec2 v1 = normalize( sp2 - sp1 );
    vec2 v2 = normalize( sp3 - sp2 );

    // determine the normal of each of the 3 segments (previous, current, next)
    vec2 n0 = vec2( -v0.y, v0.x );
    vec2 n1 = vec2( -v1.y, v1.x );
    vec2 n2 = vec2( -v2.y, v2.x );

    // determine miter lines by averaging the normals of the 2 segments
    vec2 miter_a = normalize(n0 + n1);	// miter at start of current segment
    vec2 miter_b = normalize(n1 + n2);	// miter at end of current segment

    // determine the length of the miter by projecting it onto normal and then inverse it
    float length_a = line_width / dot(miter_a, n1);
    float length_b = line_width / dot(miter_b, n1);
    
    // prevent excessively long miters at sharp corners
    // TODO: Possibly do bevels?
    if( dot(v0,v1) < -miter_limit ) {
      miter_a = n1;
      length_a = line_width;
      
      // close the gap
      if( dot(v0,n1) > 0 ) {
        gl_Position = to_clip_space( vec4( (sp1 + line_width * n0), p1.zw ) );
        g_color = v_color[0];
        g_object_id = v_object_id[0];
        EmitVertex();

        gl_Position = to_clip_space( vec4( (sp1 + line_width * n1), p1.zw ) );
        g_color = v_color[0];
        g_object_id = v_object_id[0];
        EmitVertex();
        
        gl_Position = to_clip_space( vec4( sp1, p1.zw ) );
        g_color = v_color[0];
        g_object_id = v_object_id[0];
        EmitVertex();

        EndPrimitive();
      }
      else 
      {
        gl_Position = to_clip_space( vec4( (sp1 - line_width * n1), p1.zw ) );
        g_color = v_color[0];
        g_object_id = v_object_id[0];
        EmitVertex();

        gl_Position = to_clip_space( vec4( (sp1 - line_width * n0), p1.zw ) );
        g_color = v_color[0];
        g_object_id = v_object_id[0];
        EmitVertex();

        gl_Position = to_clip_space( vec4( sp1, p1.zw ) );
        g_color = v_color[0];
        g_object_id = v_object_id[0];
        EmitVertex();

        EndPrimitive();
      }
    }

    if( dot(v1,v2) < -miter_limit ) {
      miter_b = n1;
      length_b = line_width;
    }

    // generate the triangle strip
    gl_Position = to_clip_space( vec4( (sp1 + length_a * miter_a), p1.zw ) );
    g_color = v_color[1];
    g_object_id = v_object_id[0];
    EmitVertex();

    gl_Position = to_clip_space( vec4( (sp1 - length_a * miter_a), p1.zw ) );
    g_color = v_color[1];
    g_object_id = v_object_id[0];
    EmitVertex();

    gl_Position = to_clip_space( vec4( (sp2 + length_b * miter_b), p2.zw ) );
    g_color = v_color[2];
    g_object_id = v_object_id[0];
    EmitVertex();

    gl_Position = to_clip_space( vec4( (sp2 - length_b* miter_b), p2.zw ) );
    g_color = v_color[2];
    g_object_id = v_object_id[0];
    EmitVertex();

    EndPrimitive();
  } );


char * gs_lines = (char *) SHADER_HEAD STR(
  layout( lines ) in;
  layout( triangle_strip, max_vertices = 4 ) out;

  in vec4 v_color[];
  flat in float v_object_id[];
  out vec4 g_color;
  flat out float g_object_id;

  uniform float	line_width;		// the thickness of the line in pixels
  uniform vec2	win_size;		  // the size of the viewport in pixels

  vec2 to_screen_space(vec4 pos)
  {
    // Division by w takes us to normalized device coordinates and then
    // win scale is simply the viewport transform:
    // http://www.songho.ca/opengl/gl_transform.html
    return vec2( pos.xy / pos.w ) * win_size;
  }

  vec4 to_clip_space( vec4 pos )
  {
    return vec4( (pos.xy / win_size) * pos.w , pos.zw );
  }

  void main(void)
  {
    // get the four vertices passed to the shader:
    vec4 p0 = gl_in[0].gl_Position; 
    vec4 p1 = gl_in[1].gl_Position; 

    //move to screen space
    vec2 sp0 = to_screen_space( p0 );	// start of current segment
    vec2 sp1 = to_screen_space( p1 );	// end of current segment

    // determine the direction 
    vec2 v = normalize( sp1 - sp0 );

    // determine the normal
    vec2 n = vec2( -v.y, v.x );

    // generate the triangle strip
    gl_Position = to_clip_space( vec4( (sp0 + line_width * n), p0.zw ) );
    g_color = v_color[0];
    g_object_id = v_object_id[0];
    EmitVertex();

    gl_Position = to_clip_space( vec4( (sp0 - line_width * n), p0.zw ) );
    g_color = v_color[0];
    g_object_id = v_object_id[0];
    EmitVertex();

    gl_Position = to_clip_space( vec4( (sp1 + line_width * n), p1.zw ) );
    g_color = v_color[1];
    g_object_id = v_object_id[1];
    EmitVertex();

    gl_Position = to_clip_space( vec4( (sp1 - line_width * n), p1.zw ) );
    g_color = v_color[1];
    g_object_id = v_object_id[1];
    EmitVertex();

    EndPrimitive();
  } );

char * fs_lines = (char *) SHADER_HEAD STR(
    in vec4 g_color;
    flat in float g_object_id;
    layout (location = 0) out vec4 frag_color;
    layout (location = 1) out vec4 frag_id;

    uniform bool use_uniform_color;
    uniform vec3 uniform_color;
    uniform int class_id;

    void main() {
      if (!use_uniform_color)
      {
        frag_color = g_color;
      }
      else
      {
        frag_color = vec4(uniform_color, 1.0f);
      }
      frag_id = vec4(class_id, g_object_id, 0.0, 1.0);
    });

////////////////////////////////////////////////////////////////////////////////
// Implementations
////////////////////////////////////////////////////////////////////////////////

void rendering::
init( RenderingState * render, 
      const SequenceData * seq_data,
      const Options * opts,
      const bsc::vec4i viewport )
{
  // some opengl options
  glEnable(GL_PROGRAM_POINT_SIZE);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  bsc::create_shader_prog_from_source( vs_pointcloud, 
                                       fs_pointcloud, 
                                       render->pointcloud_shader );

  bsc::create_shader_prog_from_source( vs_default, 
                                       fs_default, 
                                       render->default_shader);

  bsc::create_shader_prog_from_source( vs_fullscreen, 
                                       fs_fullscreen, 
                                       render->fullscreen_shader);

  bsc::create_shader_prog_from_source( vs_lines,
                                       gs_line_adj,
                                       fs_lines,
                                       render->line_adj_shader );

  bsc::create_shader_prog_from_source( vs_lines,
                                       gs_lines,
                                       fs_lines,
                                       render->line_shader );

  // initialize common stuff
  rendering::init_axes( &(render->axes_geo));

  rendering::init_frustum( &(render->frustum_geo), 
                           &(seq_data->depth_intrinsics));

  rendering::init_grid( &render->grid_geo, 100);

  rendering::init_disc( &render->disc_geo );


  rendering::init_trajectory( &(render->trajectory_geo), 
                              &(seq_data->current_xforms[0]), 
                              &(seq_data->pairwise_xforms_inliers[0]),
                              opts,
                              0,
                              opts->frame_end - 1 ); 
  render->trajectory_dirty = false;

  // if we ground truth correspondences let's initialize them too!
  if ( seq_data->n_gt_corrs > 0 )
  {
    rendering::init_gt_corrs( &(render->gt_corrs_geo),
                              seq_data->gt_corrs,
                              seq_data->n_gt_corrs,
                              &(seq_data->current_xforms[0]) );
    render->gt_corrs_dirty = false;
  }

  // initialize bboxes -> This is very similar to the pointcloud initialization.
  render->bboxes_geo= (bsc::gpu_geometry*)malloc( sizeof( bsc::gpu_geometry ) * 
                                                            opts->frame_end );
  render->bboxes_dirty = (bool*)malloc( sizeof(bool) * opts->frame_end );
  for ( i32 i = 0 ; 
        i < opts->frame_end ; 
        i++  )
  {
    render->bboxes_geo[i].vao = -1;
    render->bboxes_dirty[i] = true;
  }

  // initialize fake fullscreen geo
  glGenVertexArrays(1, &render->fullscreen_vao);
  render->viewport = bsc::vec4i(-1,-1,-1,-1);
  render->pixel_ratio = -1;
}

void rendering::
draw( RenderingState * render, 
      const SequenceData * seq_data,
      const ViewOptions * view_opts,
      const OptimizationParameters * optim_params,
      const Options * opts,
      const bsc::vec4i viewport )
{

  if( viewport != render->viewport ||
      bsc::g_window.pixel_ratio != render->pixel_ratio )
  {
    rendering::init_framebuffer( render, viewport );
  }
  // check if any of the dirty flags is set
  // TODO: Possibly streaming the drawing / mapping the buffer is better here.
  if ( render->trajectory_dirty )
  {
    rendering::init_trajectory( &(render->trajectory_geo), 
                                &(view_opts->xforms_ptr[0]),
                                &(seq_data->pairwise_xforms_inliers[0]),
                                opts,
                                view_opts->frame_idx_A,
                                view_opts->frame_idx_B );
    render->trajectory_dirty = false;
  }
  if ( seq_data->n_gt_corrs > 0 && 
       render->gt_corrs_dirty )
  {
    rendering::init_gt_corrs( &(render->gt_corrs_geo),
                                seq_data->gt_corrs,
                                seq_data->n_gt_corrs,
                                &(view_opts->xforms_ptr[0]),
                                view_opts->frame_idx_A,
                                view_opts->frame_idx_B );
    render->gt_corrs_dirty = false;
  }

  if (seq_data->reconstruction )
  {
    if ( !seq_data->cp_corrs.empty() &&
         render->cp_corrs_dirty )
    {
      rendering::init_cp_corrs( &(render->cp_corrs_geo),
                                  seq_data,
                                  view_opts );
      render->cp_corrs_dirty = false;
    }
    if ( !seq_data->geometrical_corrs.empty() &&
         seq_data->model->structure && 
         render->geometrical_corrs_dirty )
    {
      rendering::init_geometrical_corrs( &(render->geometrical_corrs_geo),
                                        seq_data,
                                        view_opts );
      render->geometrical_corrs_dirty = false;
    }
    if ( !seq_data->hierarchical_corrs.empty() &&
         seq_data->model->structure && 
         render->hierarchical_corrs_dirty )
    {
      rendering::init_hierarchical_corrs( &(render->hierarchical_corrs_geo),
                                      seq_data,
                                      view_opts );
      render->hierarchical_corrs_dirty = false;
    }

    // if all bboxes are marked clean, this just passes through
    init_bboxes( render->bboxes_geo,
                 render->bboxes_dirty,
                 seq_data,
                 opts );
  
  }

  // Start Drawing!

  // Bind frambebuffer
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, render->fullscreen_framebuffer);
  
  // Clear Window
  bsc::window & g_window = bsc::get_window();
  bsc::clear_window( viewport );
  
  bsc::mat4 vp = render->proj * render->view;
  bsc::use_program(render->line_shader);
  bsc::set_uniform(render->line_shader, "line_width", view_opts->line_width );
  bsc::set_uniform(render->line_shader, "win_size", bsc::vec2( 1024, 768 ) );
  bsc::set_uniform(render->line_shader, "use_uniform_color", false);
  bsc::set_uniform(render->line_shader, "mvp", vp);

  if ( seq_data->reconstruction && 
       !seq_data->cp_corrs.empty() && 
       view_opts->show_cp_corrs )
  {
    bsc::set_uniform(render->line_shader, "class_id", CP_CORR_ID);
    draw( &render->cp_corrs_geo, GL_LINES);
  }

  if ( !seq_data->geometrical_corrs.empty() && 
       view_opts->show_geometrical_corrs )
  {
    bsc::set_uniform(render->line_shader, "class_id", GEOMETRICAL_CORR_ID);
    draw( &render->geometrical_corrs_geo, GL_LINES);
  }

  if ( !seq_data->hierarchical_corrs.empty() && 
       view_opts->show_hierarchical_corrs )
  {
    bsc::set_uniform(render->line_shader, "class_id", HIERARCHICAL_CORR_ID);
    draw( &render->hierarchical_corrs_geo, GL_LINES);
  }

  bsc::use_program(render->line_adj_shader);
  bsc::set_uniform(render->line_adj_shader, "line_width", view_opts->line_width );
  bsc::set_uniform(render->line_adj_shader, "miter_limit", 0.75f );
  bsc::set_uniform(render->line_adj_shader, "win_size", bsc::vec2( 1024, 768 ) );
  bsc::set_uniform(render->line_adj_shader, "use_uniform_color", false);
  bsc::set_uniform(render->line_adj_shader, "mvp", vp);

  if (view_opts->show_axes)
  {
    bsc::set_uniform(render->line_adj_shader, "class_id", 0);
    draw( &render->axes_geo, GL_LINES_ADJACENCY );
  }

  if (view_opts->show_trajectory)
  {
    bsc::set_uniform(render->line_adj_shader, "class_id", TRAJECTORY_ID);
    draw( &render->trajectory_geo, GL_LINES_ADJACENCY );
  }

  bsc::use_program(render->default_shader);
  bsc::set_uniform(render->default_shader, "use_uniform_color", false);
  bsc::set_uniform(render->default_shader, "mvp", vp);
  bsc::set_uniform(render->default_shader, "class_id", 0);

  bsc::set_uniform(render->default_shader, "use_uniform_color", true);
  if (view_opts->show_grid)
  {
    bsc::set_uniform(render->default_shader, "uniform_color", 
                                              bsc::vec3( 0.6f, 0.6f, 0.6f ) );
    draw( &render->grid_geo, GL_LINES );
  }


  for (i32 frame_idx = 0;
       frame_idx < opts->frame_end;
       frame_idx++ )
  {
    bool should_skip = view_opts->show_only_pair ? ( frame_idx != view_opts->frame_idx_A && frame_idx != view_opts->frame_idx_B )
                                                 : ( frame_idx < view_opts->frame_idx_A  || frame_idx > view_opts->frame_idx_B );
    if ( view_opts->frame_skip != 0 && 
         frame_idx % view_opts->frame_skip != 0 ) continue;
    if ( should_skip ) continue;

    // Convert from double to float matrix... how to do implicit conversions?
    bsc::mat4d tmp = view_opts->xforms_ptr[ frame_idx ];
    bsc::mat4 model;
    model[0] = bsc::vec4( tmp[0][0], tmp[0][1], tmp[0][2], tmp[0][3] );
    model[1] = bsc::vec4( tmp[1][0], tmp[1][1], tmp[1][2], tmp[1][3] );
    model[2] = bsc::vec4( tmp[2][0], tmp[2][1], tmp[2][2], tmp[2][3] );
    model[3] = bsc::vec4( tmp[3][0], tmp[3][1], tmp[3][2], tmp[3][3] );


    bsc::vec3i col = qual_map[frame_idx % 9];
    bsc::use_program(render->default_shader);

    bsc::set_uniform(render->default_shader, "use_uniform_color", true);
    bsc::set_uniform(render->default_shader, "uniform_color", 
                                                   bsc::vec3( col.r / 255.0f,
                                                              col.g / 255.0f,
                                                              col.b / 255.0f));
    bsc::set_uniform(render->default_shader, "class_id", 0);
    if (view_opts->show_cameras)
    {
      bsc::set_uniform(render->default_shader, "mvp", 
                                          render->proj * 
                                          render->view * model);
    
      bsc::draw(&render->frustum_geo, GL_LINE_STRIP);
    }

    if ( seq_data->reconstruction && view_opts->show_bboxes)
    {
      bsc::set_uniform(render->default_shader, "mvp", 
                                                  render->proj * render->view );
      bsc::draw(&(render->bboxes_geo[frame_idx]), GL_LINES );
    }

    if (view_opts->show_frames && 
        render->pointclouds_geo[ frame_idx ].vao != -1 )
    {
      // TODO: Test speed of doing mvp in shader, versus passing it
      bsc::use_program(render->pointcloud_shader);
      bsc::set_uniform(render->pointcloud_shader, "max_dist", 
                                                          view_opts->max_dist );
      bsc::set_uniform(render->pointcloud_shader, "max_height", 
                                                        view_opts->max_height );
      bsc::set_uniform(render->pointcloud_shader, "point_size", 
                       (r32)bsc::g_window.pixel_ratio * view_opts->point_size );
      bsc::set_uniform(render->pointcloud_shader, "coloring_mode", 
                                               (int)view_opts->feat_disp_mode );
      bsc::set_uniform(render->pointcloud_shader, "mv", 
                         render->view * model );
      bsc::set_uniform(render->pointcloud_shader, "m",  
                                        model );
      bsc::set_uniform(render->pointcloud_shader, "p", render->proj);
      bsc::set_uniform(render->pointcloud_shader, "class_id", POINTCLOUD_ID);
      bsc::set_uniform(render->pointcloud_shader, "object_id", frame_idx );
      draw( &(render->pointclouds_geo[ frame_idx ]), GL_POINTS);
    }
  }

  bsc::use_program(render->default_shader);
  bsc::set_uniform(render->default_shader, "use_uniform_color", true );
  bsc::set_uniform(render->default_shader, "uniform_color", 
                                                   bsc::vec3( 1.0f,
                                                              1.0f,
                                                              1.0f ) );
  bsc::set_uniform(render->default_shader, "class_id", PROXY_ID);
  
  if ( view_opts->show_proxies && 
       seq_data->model && 
       seq_data->model->structure )
  {
    for ( i32 proxy_idx = 0 ;
          proxy_idx < seq_data->model->structure->NProxies() ;
          proxy_idx++ )
    {
      if( view_opts->inspection_mode && 
          view_opts->selected_class_id == PROXY_ID &&
          view_opts->selected_object_id != proxy_idx ) continue;
          
      Proxy * proxy = seq_data->model->structure->GetProxy( proxy_idx );
     
      if ( view_opts->proxy_level > -1 && 
           view_opts->proxy_level != proxy->level )
      {
        continue;
      }
      std::vector<i32> shapes_ids;
      optproxy::get_shape_ind_connected_to_proxy( proxy, shapes_ids );
      i32 should_show = 0; 
      for ( i32 idx = 0;
            idx < shapes_ids.size();
            idx++ )
      {
        if ( shapes_ids[ idx ] <= view_opts->frame_idx_B && 
             shapes_ids[ idx ] >= view_opts->frame_idx_A )
        {
          should_show = 1;
          break;
        }
      }
      if ( !should_show ) 
      {
        continue;
      }

      // check if any shape contains the proxy.
      FETFeature * feature = proxy->feature;
      FETShape * shape = feature->Shape();
      R3Affine T = shape->Transformation( view_opts->xform_scheme );
      R3Point pos = feature->Position();
      R3Vector normal = feature->Normal();
      pos.Transform( T );
      normal.Transform( T );

      bsc::vec3 p( pos.X(), pos.Y(), pos.Z() );
      bsc::vec3 n( normal.X(), normal.Y(), normal.Z() );
      p += (proxy->level+1) * view_opts->proxy_offset * (-n);

      bsc::vec3 v = bsc::normalize( bsc::cross( n, bsc::vec3( 0.0f, 1.0f, 0.0f ) ) );
      bsc::vec3 u = bsc::normalize( bsc::cross( v, n ) );

      bsc::mat4 xform;
      xform[0] = bsc::vec4( u, 0.0f );
      xform[1] = bsc::vec4( n, 0.0f );
      xform[2] = bsc::vec4( v, 0.0f );
      xform[3] = bsc::vec4( p, 1.0f );
      
      bsc::vec3i col;
      if ( view_opts->proxy_disp_mode == SHOW_LEVEL )
        col = qual_map[proxy->level % 9];
      if ( view_opts->proxy_disp_mode == SHOW_PROXY_IDX )
        col = qual_map[proxy_idx % 9];
      if ( view_opts->proxy_disp_mode == SHOW_PARENT_IDX )
      {
        i32 structure_idx = proxy->structure_index;
        if ( proxy->level == 0 )
        {
          Proxy * parent = proxy->parent;
          if ( parent != NULL )
            structure_idx = parent->structure_index;
        }
        col = divergent_map[structure_idx % 15];
      }
      bsc::set_uniform( render->default_shader, "object_id", proxy_idx );
      bsc::set_uniform(render->default_shader, "uniform_color", 
                                                   bsc::vec3( col.r / 255.0f,
                                                              col.g / 255.0f,
                                                              col.b / 255.0f));
      bsc::set_uniform(render->default_shader, "mvp",  vp * xform );

      bsc::draw( &(render->disc_geo), GL_TRIANGLE_FAN );
    }
  }

  vp = render->proj * render->view;
  bsc::use_program(render->line_shader);
  bsc::set_uniform(render->line_shader, "line_width", view_opts->line_width );
  bsc::set_uniform(render->line_shader, "win_size", bsc::vec2( 1024, 768 ) );
  bsc::set_uniform(render->line_shader, "use_uniform_color", false);
  bsc::set_uniform(render->line_shader, "mvp", vp);

  if (seq_data->n_gt_corrs > 0 && 
      view_opts->show_gt_corrs)
  {
    glDisable( GL_DEPTH_TEST );
    draw(&render->gt_corrs_geo, GL_LINES);
    glEnable (GL_DEPTH_TEST );
  }
  

// RENDER TO SCREEN!
  glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
  
  glDisable( GL_DEPTH_TEST );
  glViewport( g_window.pixel_ratio * viewport[0], 
              g_window.pixel_ratio * viewport[1], 
              g_window.pixel_ratio * viewport[2], 
              g_window.pixel_ratio * viewport[3] );
  glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
  
  bsc::use_program( render->fullscreen_shader );
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, render->color_buffer );
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, render->index_buffer );
  bsc::set_uniform(render->fullscreen_shader, "u_color", 0);
  bsc::set_uniform(render->fullscreen_shader, "u_index", 1);


  glBindVertexArray( render->fullscreen_vao );
  glDrawArrays(GL_TRIANGLES, 0, 3);
  glBindVertexArray(0);

  glEnable( GL_DEPTH_TEST );

}

void rendering::
request_update( RenderingState * render, 
                const Options * opts )
{
    render->trajectory_dirty           = true;
    render->gt_corrs_dirty             = true;
    render->cp_corrs_dirty             = true;
    render->geometrical_corrs_dirty       = true;
    render->hierarchical_corrs_dirty         = true;
    
    for ( i32 i = 0; i < opts->frame_end; i++ )
    {
      render->bboxes_dirty[i] = true;
    }
}

void rendering::
init_framebuffer( RenderingState *render, bsc::vec4i viewport )
{
  // modify viewport
  int w = bsc::g_window.pixel_ratio * (viewport.z - viewport.x);
  int h = bsc::g_window.pixel_ratio * (viewport.w - viewport.y);
  render->viewport = viewport;
  
  render->pixel_ratio = bsc::g_window.pixel_ratio;

  if(!render->fullscreen_framebuffer) { glDeleteFramebuffers(1, &render->fullscreen_framebuffer); }
  if(!render->color_buffer) { glDeleteTextures(1, &render->color_buffer); }
  if(!render->index_buffer) { glDeleteTextures(1, &render->index_buffer); }
  if(!render->depth_renderbuffer) { glDeleteRenderbuffers(1, &render->depth_renderbuffer); }

  glGenFramebuffers(1, &render->fullscreen_framebuffer );
  glGenTextures( 1, &(render->color_buffer) );
  glGenTextures( 1, &(render->index_buffer ) );
  glGenRenderbuffers( 1, &render->depth_renderbuffer );

  glBindFramebuffer( GL_FRAMEBUFFER, render->fullscreen_framebuffer );

  // Color
  glBindTexture( GL_TEXTURE_2D, render->color_buffer );
  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  glFramebufferTexture( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, render->color_buffer, 0 );

  // Index
  glBindTexture( GL_TEXTURE_2D, render->index_buffer );
  glTexImage2D( GL_TEXTURE_2D, 0, GL_RG32F, w, h,  0, GL_RG, GL_FLOAT, NULL );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  glFramebufferTexture( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, render->index_buffer, 0 );

  glReadBuffer(GL_NONE);
  GLuint attachments[2] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
  glDrawBuffers(2, attachments);

  // Depth
  glBindRenderbuffer(GL_RENDERBUFFER, render->depth_renderbuffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w, h );
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, render->depth_renderbuffer);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
  {
    printf("Framebuffer not complete!\n");
  }
  glBindFramebuffer( GL_FRAMEBUFFER, 0 );
}

i32 rendering::
init_pointclouds_from_frames( RenderingState * render,
                              const SequenceData * seq_data,
                              const Options * opts )
{
  i32 n_points = 0;
  
  static bsc::img color_img;
  static bsc::img_u16 depth_img;

  // memory init
  render->pointclouds_geo = (bsc::gpu_geometry*)malloc( sizeof( bsc::gpu_geometry ) * 
                                                        opts->frame_end );
  render->pointclouds_dirty = (bool*)malloc( sizeof(bool) * 
                                             opts->frame_end );
  print_mat3( seq_data->depth_intrinsics );
  for ( i32 i = 0 ; i < opts->frame_end ; i++ )
  {
    render->pointclouds_geo[i].vao = -1;
    render->pointclouds_dirty[i] = true;
  }

  // read frames
  bsc::stop_watch timer;
  timer.start();
  if (opts->print_verbose)
  {
    printf("Starting to read the frames...\n");
  }

  for (i32 frame_idx = 0;
       frame_idx < opts->frame_end;
       frame_idx++ )
  {
    // Report status
    if ( opts->print_verbose )
    {
      if (frame_idx > 0) { printf("\r"); }
      printf("Frame %d/%d", frame_idx + 1, opts->frame_end);
      fflush(stdout);
    }

    // Prepare names
    char color_filename[512], depth_filename[512];
    sprintf(color_filename, "%s/%s", seq_data->color_dir,
            seq_data->color_names[frame_idx]);
    sprintf(depth_filename, "%s/%s", seq_data->depth_dir,
            seq_data->depth_names[frame_idx]);

    // Read in images
    color_img.read( color_filename );
    depth_img.read( depth_filename );

    bool bitshift = !strcmp(seq_data->dataset, "sun3d");

    n_points += rendering::init_pointcloud(
                                          &(render->pointclouds_geo[frame_idx]),
                                          &color_img,
                                          &depth_img,
                                          seq_data->depth_intrinsics,
                                          opts->downsample_factor,
                                          bitshift );
    
    render->pointclouds_dirty[ frame_idx ] = false;
  }

  if ( opts->print_verbose )
  {
    printf( "\nTook %f sec.\n", timer.elapsed() );
  }
  return n_points;
}

i32 rendering::
init_pointclouds_from_reconstruction( RenderingState * render,
                                      const SequenceData * seq_data,
                                      const ViewOptions * view_opts,
                                      const Options * opts )
{
  // useful variables
  const FETReconstruction * recstr = seq_data->reconstruction;
  i32 n_points = 0;

  // gpu stroage init
  render->pointclouds_geo = (bsc::gpu_geometry*)malloc( sizeof( bsc::gpu_geometry ) * 
                                                  opts->frame_end );
  render->pointclouds_dirty = (bool*)malloc( sizeof(bool) * 
                                             opts->frame_end );
  for ( i32 i = 0 ; i < opts->frame_end ; i++ )
  {
    render->pointclouds_geo[i].vao = -1;
    render->pointclouds_dirty[i] = true;
  }


  // temp local storage
  bsc::geometry_data pointcloud;
  i32 base_size = seq_data->depth_res.x * seq_data->depth_res.y;
  pointcloud.positions = (r32*)malloc( sizeof(bsc::vec3) * 3 * base_size );
  pointcloud.normals = (r32*)malloc( sizeof(bsc::vec3) * 3 * base_size );
  pointcloud.colors_a = (u8*)malloc( sizeof(u8) * 4 * base_size );
  pointcloud.colors_b = (u8*)malloc( sizeof(u8) * 4 * base_size );
  
  // copy data!
  for ( i32 shape_idx = 0 ;
        shape_idx < opts->frame_end ;
        shape_idx++ )
  {
    FETShape * shape = recstr->Shape( shape_idx );
    pointcloud.n_vertices = shape->NFeatures(); 

    if ( shape->NFeatures() == 0 ) continue;

    i64 pos_idx = 0, nor_idx = 0, col_a_idx = 0, col_b_idx = 0;
    for ( i32 feat_idx = 0 ;
          feat_idx < shape->NFeatures() ;
          feat_idx++ )
    {
      FETFeature * f = shape->Feature( feat_idx );
      R3Point p = f->Position();
      R3Vector n = f->Normal();
      RNRgb c = f->Color();

      pointcloud.positions[ pos_idx++ ] = p.X();
      pointcloud.positions[ pos_idx++ ] = p.Y();
      pointcloud.positions[ pos_idx++ ] = p.Z();

      pointcloud.normals[ nor_idx++ ] = n.X();
      pointcloud.normals[ nor_idx++ ] = n.Y();
      pointcloud.normals[ nor_idx++ ] = n.Z();

      pointcloud.colors_a[ col_a_idx++ ] = (u8)(c.R() * 255.0f);
      pointcloud.colors_a[ col_a_idx++ ] = (u8)(c.G() * 255.0f);
      pointcloud.colors_a[ col_a_idx++ ] = (u8)(c.B() * 255.0f);
      pointcloud.colors_a[ col_a_idx++ ] = f->GeneratorType(); // For coloring!

      // Issue - the marker will change.
      i32 marker = f->primitive_marker;
      pointcloud.colors_b[ col_b_idx++ ] = (u8)qual_map[ marker % 9 ].r ;
      pointcloud.colors_b[ col_b_idx++ ] = (u8)qual_map[ marker % 9 ].g;
      pointcloud.colors_b[ col_b_idx++ ] = (u8)qual_map[ marker % 9 ].b;
      pointcloud.colors_b[ col_b_idx++ ] = 255;

      n_points++;
    }
    bsc::init_gpu_geo( &pointcloud, &(render->pointclouds_geo[shape_idx]), 
                       POSITION | NORMAL | COLOR_A | COLOR_B );
    render->pointclouds_dirty[ shape_idx ] = false;
  }

  free( pointcloud.positions );
  free( pointcloud.normals );
  free( pointcloud.colors_a );
  free( pointcloud.colors_b );

  return n_points;
}

void rendering::
init_axes( bsc::gpu_geometry *cmd )
{
  bsc::geometry_data axes;
  r32 positions[36] = {-1,  0,  0,
                        0,  0,  0,
                        3,  0,  0,
                        4,  0,  0,

                        0, -1,  0,
                        0,  0,  0,
                        0,  3,  0,
                        0,  4,  0,
                        
                        0,  0, -1,
                        0,  0,  0,
                        0,  0,  3,
                        0,  0,  4 };

  u8 colors[48] = {255, 0, 0, 255, 255, 0, 0, 255,
                   255, 0, 0, 255, 255, 0, 0, 255,
                   0, 255, 0, 255, 0, 255, 0, 255,
                   0, 255, 0, 255, 0, 255, 0, 255,
                   0, 0, 255, 255, 0, 0, 255, 255,
                   0, 0, 255, 255, 0, 0, 255, 255};

  u32 indices[12] = { 0, 1, 2, 3, 
                      4, 5, 6, 7, 
                      8, 9, 10, 11 };
  
  axes.positions = &(positions[0]);
  axes.colors_a  = &(colors[0]);
  axes.indices =   &(indices[0]);
  axes.n_vertices = 12;
  axes.n_elements = 12;

  bsc::init_gpu_geo( &axes, cmd, POSITION | COLOR_A | STRUCTURED );
}

void rendering::
init_bboxes( bsc::gpu_geometry *bboxes_geo,
             bool *dirty,
             const SequenceData *seq_data,
             const Options *opts )
{

  for ( i32 bbox_idx = 0;
        bbox_idx < opts->frame_end;
        bbox_idx++  )
  {
    if ( !dirty[ bbox_idx ] )
    {
      continue;
    }

    // TODO: This is a lot of gpu memory moving. Should be a better way
    if ( bboxes_geo[bbox_idx].vao != -1 )
    {
      bsc::free_gpu_geo( &(bboxes_geo[bbox_idx] ) );
    }

    const R3Box box = seq_data->bboxes[bbox_idx];//recstr->Shape(bbox_idx)->BBox();

    // extract 8 vertices
    bsc::vec3 p_a( box.Min().X(), box.Min().Y(), box.Min().Z() );
    bsc::vec3 p_b( box.Max().X(), box.Min().Y(), box.Min().Z() );
    bsc::vec3 p_c( box.Max().X(), box.Min().Y(), box.Max().Z() );
    bsc::vec3 p_d( box.Min().X(), box.Min().Y(), box.Max().Z() );
    
    bsc::vec3 p_e( box.Min().X(), box.Max().Y(), box.Min().Z() );
    bsc::vec3 p_f( box.Max().X(), box.Max().Y(), box.Min().Z() );
    bsc::vec3 p_g( box.Max().X(), box.Max().Y(), box.Max().Z() );
    bsc::vec3 p_h( box.Min().X(), box.Max().Y(), box.Max().Z() );
    
    //put them in an array as lines
    bsc::vec3 positions[24] = 
    { 
      // bottom edges
      p_a, p_b,
      p_b, p_c,
      p_c, p_d,
      p_d, p_a,

      // top edges
      p_e, p_f,
      p_f, p_g,
      p_g, p_h,
      p_h, p_e,

      // side edges
      p_a, p_e,
      p_b, p_f,
      p_c, p_g,
      p_d, p_h
    };

    bsc::geometry_data bbox_data;
    bbox_data.positions = &(positions[0].x);
    bbox_data.n_vertices = 24;

    bsc::init_gpu_geo( &bbox_data, &(bboxes_geo[bbox_idx]), POSITION ); 
    dirty[bbox_idx] = false;
  }

}

void rendering::
init_trajectory( bsc::gpu_geometry *traj_geo,
                 const bsc::mat4d *extrinsics,
                 const int *n_inliers,
                 const Options *opts,
                 const i32 low_lim, 
                 const i32 high_lim )
{
  using namespace bsc;
  // NOTE: This is inefficient, but does not happen too often
  if ( traj_geo->vao != -1 )
  {
    free_gpu_geo( traj_geo );
  }

  i32 n = (high_lim - low_lim) + 1;
  if ( n <= 0 ) return;
  r32 positions[ 3 * (n + 2) ];
  u8 colors[ 4 * (n + 2) ];
  u32 indices[ 4 * n ];

  vec3d start = 2.0 * ( vec3d(extrinsics[0][3])   - vec3d(extrinsics[1][3]) );
  vec3d end  = 2.0 * ( vec3d(extrinsics[n-1][3]) - vec3d(extrinsics[n-2][3]) );

  // start / end
  positions[0] = start.x; positions[1] = start.y; positions[2] = start.z;
  positions[3*n+3] = end.x; positions[3*n+4] = end.y; positions[3*n+5] = end.z;
  colors[0] = 255; colors[1] = 255; colors[2] = 255; colors[3] = 255;
  colors[4*n+4] = 255; colors[4*n+5] = 255; colors[4*n+6] = 255; colors[4*n+7] = 255;

  // store positions
  i32 loc_idx = 3, col_idx = 4;
  i32 n_pos = 0;
  for ( i32 i = low_lim ; i <= high_lim ; i++  )
  {
    vec3d pos = vec3d( extrinsics[i][3] );
    positions[loc_idx++] = pos.x;
    positions[loc_idx++] = pos.y;
    positions[loc_idx++] = pos.z;

    r32 t = (r32)i / (r32)opts->frame_end;
    i32 col_idx_A = floor(8 * t);
    i32 col_idx_B = std::min( col_idx_A + 1, 8 );
    float weight = mmath::sigmoid( n_inliers[i], 5.0f, 0.5f );
    bsc::vec3 col_A( transition_map[ col_idx_A ].r, 
                     transition_map[ col_idx_A ].g, 
                     transition_map[ col_idx_A ].b );
    bsc::vec3 col_B( transition_map[ col_idx_B ].r, 
                     transition_map[ col_idx_B ].g, 
                     transition_map[ col_idx_B ].b );
    r32 t2 = 1.0f - ((8.0f * t) - col_idx_A);

    bsc::vec3 col = weight * (t2 * col_A + (1.0f - t2) * col_B);

    colors[ col_idx++ ] = col.r;
    colors[ col_idx++ ] = col.g;
    colors[ col_idx++ ] = col.b;
    colors[ col_idx++ ] = 255;
    n_pos++;
  }

  i32 ind_idx = 0;
  for ( int i = 1 ; i < n + 1 ; ++i )
  {
    indices[ ind_idx++ ] = i - 1;
    indices[ ind_idx++ ] = i;
    indices[ ind_idx++ ] = i + 1;
    indices[ ind_idx++ ] = i + 2;
  }
  
  bsc::geometry_data traj_data;
  traj_data.positions  = &(positions[0]);
  traj_data.colors_a   = &(colors[0]);
  traj_data.indices    = &(indices[0]);
  traj_data.n_vertices = n + 2;
  traj_data.n_elements = 4*n;

  bsc::init_gpu_geo( &traj_data, traj_geo, POSITION | COLOR_A | STRUCTURED ); 
}

void rendering::
init_gt_corrs( bsc::gpu_geometry * geo,
               const Correspondence * correspondences,
               const i32 n_corrs,
               const bsc::mat4d * xforms,
               const i32 lim_low, const i32 lim_high )
{
  if ( geo->vao != -1 ) 
  {
    bsc::free_gpu_geo( geo );
  }

  bsc::geometry_data corrs;
  corrs.positions  = (r32*)malloc( 6 * n_corrs * sizeof(r32)  );
  corrs.colors_a   = (u8*)malloc( 8 * n_corrs * sizeof(u8) ); 
  corrs.n_vertices = 2 * n_corrs;
  
  i32 min_dist = 1e9;
  i32 max_dist = -1e9;
  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        corr_idx++ )
  {
    i32 idx_A = correspondences[corr_idx].idx_A;
    i32 idx_B = correspondences[corr_idx].idx_B;
    i32 diff = abs( idx_A - idx_B );
    if ( diff > max_dist ) max_dist = diff;
    if ( diff < min_dist ) min_dist = diff;
  }

  i32 pos_idx = 0, col_idx = 0;
  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        corr_idx++ )
  {
    i32 idx_A = correspondences[corr_idx].idx_A;
    i32 idx_B = correspondences[corr_idx].idx_B;
    i32 diff = abs( idx_A - idx_B );


    if ( idx_A < lim_low || idx_A > lim_high ||
         idx_B < lim_low || idx_B > lim_high )
    {
      continue;
    }

    bsc::vec4d pos_A  = bsc::vec4d( correspondences[corr_idx].point_A, 1.0f );
    bsc::vec4d pos_B  = bsc::vec4d( correspondences[corr_idx].point_B, 1.0f );

    pos_A = xforms[idx_A] * pos_A; 
    pos_B = xforms[idx_B] * pos_B; 

    r32 t = (r32)(diff - min_dist) / (r32)(max_dist - min_dist);
    i32 col_idx_A = round(t);
    i32 col_idx_B = col_idx_A + 1;
    bsc::vec3 col_A( red2blue[ col_idx_A ].r, 
                     red2blue[ col_idx_A ].g, 
                     red2blue[ col_idx_A ].b );
    bsc::vec3 col_B( red2blue[ col_idx_B ].r, 
                     red2blue[ col_idx_B ].g, 
                     red2blue[ col_idx_B ].b );
    r32 t2 = ((2.0f * t) - col_idx_A);
    bsc::vec3 col = (1.0f - t2) * col_A + t2 * col_B; 

    // copy in data
    corrs.positions[pos_idx++] = pos_A.x; 
    corrs.positions[pos_idx++] = pos_A.y; 
    corrs.positions[pos_idx++] = pos_A.z; 
                                           
    corrs.positions[pos_idx++] = pos_B.x; 
    corrs.positions[pos_idx++] = pos_B.y; 
    corrs.positions[pos_idx++] = pos_B.z; 
                                           
    corrs.colors_a[col_idx++]  = col.r;
    corrs.colors_a[col_idx++]  = col.g;
    corrs.colors_a[col_idx++]  = col.b;
    corrs.colors_a[col_idx++]  = 255;

    corrs.colors_a[col_idx++]  = col.r;
    corrs.colors_a[col_idx++]  = col.g;
    corrs.colors_a[col_idx++]  = col.b;
    corrs.colors_a[col_idx++]  = 255;
  }

  bsc::init_gpu_geo(&corrs, geo, POSITION | COLOR_A);

  free( corrs.positions );
  free( corrs.colors_a );

}

void rendering::
init_geometrical_corrs( bsc::gpu_geometry * geo,
                      const SequenceData * seq_data,
                      const ViewOptions * view_opts )
{
  if ( geo->vao != -1 ) 
  {
    bsc::free_gpu_geo( geo );
  }

  const std::unordered_map<FETShape*, Proxy*> * shape_to_proxy = 
                                                seq_data->model->shape_to_proxy;
  const std::vector<FETCorrespondence*> & corrs = seq_data->geometrical_corrs;
  i32 n_corrs = corrs.size(); 

  bsc::geometry_data corrs_geo;
  corrs_geo.positions  = (r32*)malloc( 6 * n_corrs * sizeof(r32)  );
  corrs_geo.colors_a   = (u8*)malloc( 8 * n_corrs * sizeof(u8) ); 
  corrs_geo.n_vertices = 0;
  i32 pos_idx = 0, col_idx = 0;

  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        corr_idx ++ )
  {
    FETCorrespondence * corr = corrs[ corr_idx ];
    if ( !corr ) { printf("To Model corr invalid\n"); continue; }
    FETFeature * f1 = corr->Feature(0);
    FETFeature * f2 = corr->Feature(1);
    if ( !f1 || !f2 ) { printf("To Model corr features invalid\n"); continue; }
    FETShape * s1 = f1->Shape();
    FETShape * s2 = f2->Shape();
    if ( !s1 || !s2 ) { printf("To Model corr shapes invalid\n"); continue; }

    R3Point p1 = f1->Position();
    p1.Transform( s1->Transformation( view_opts->xform_scheme ) );
    R3Point p2 = f2->Position();
    p2.Transform( s2->Transformation( view_opts->xform_scheme ) );

    R3Vector n1 = f1->Normal();
    n1.Transform( s1->Transformation( view_opts->xform_scheme ) );
    n1 = -n1;
    R3Vector n2 = f2->Normal();
    n2.Transform( s2->Transformation( view_opts->xform_scheme ) );
    n2 = -n2;

    i32 l1 = 0;
    i32 l2 = 0;
    
    // bit of C++11 abuse...
    auto shape_proxy1 = shape_to_proxy->find( f1->Shape() );
    auto shape_proxy2 = shape_to_proxy->find( f2->Shape() );

    if ( shape_proxy1 != shape_to_proxy->end() ) 
    {
      l1 = shape_proxy1->second->level + 1;
    }

    if ( shape_proxy2 != shape_to_proxy->end() ) 
    {
      l2 = shape_proxy2->second->level + 1;
    }

    // if level is 0, then just check shape.
    // else get inlier shape ids, and check if they are contained.
    bool should_show = 0;
    if ( l1 == 0 && is_within( f1->Shape()->reconstruction_index, view_opts ) )
    {
      should_show |= 1;
    }
    if ( l2 == 0 && is_within( f2->Shape()->reconstruction_index, view_opts ) )
    {
      should_show |= 1;
    }
    if ( l1 > 0 && is_within( shape_proxy1->second, view_opts ) )
    {
      should_show |= 1;
    }
    if ( l2 > 0 && is_within( shape_proxy2->second, view_opts ) )
    {
      should_show |= 1;
    }
    if ( !should_show ) continue;


    // Note: How to figure out level?!
    p1 += l1 * view_opts->proxy_offset * n1;
    p2 += l2 * view_opts->proxy_offset * n2;

    bsc::vec3i col_A = qual_map[ std::max( l1-1, l2-1) % 9 ];
    bsc::vec3i col_B = qual_map[ std::max( l1-1, l2-1) % 9 ];

    // copy in data
    corrs_geo.positions[pos_idx++] = p1.X(); 
    corrs_geo.positions[pos_idx++] = p1.Y(); 
    corrs_geo.positions[pos_idx++] = p1.Z(); 
                                          
    corrs_geo.positions[pos_idx++] = p2.X(); 
    corrs_geo.positions[pos_idx++] = p2.Y(); 
    corrs_geo.positions[pos_idx++] = p2.Z(); 
                                          
    corrs_geo.colors_a[col_idx++]  = col_A.r;
    corrs_geo.colors_a[col_idx++]  = col_A.g;
    corrs_geo.colors_a[col_idx++]  = col_A.b;
    corrs_geo.colors_a[col_idx++]  = 255;

    corrs_geo.colors_a[col_idx++]  = col_B.r;
    corrs_geo.colors_a[col_idx++]  = col_B.g;
    corrs_geo.colors_a[col_idx++]  = col_B.b;
    corrs_geo.colors_a[col_idx++]  = 255;

    corrs_geo.n_vertices += 2;
  }

  if ( corrs_geo.n_vertices > 0 )
  {
    bsc::init_gpu_geo(&corrs_geo, geo, POSITION | COLOR_A);
  }

  free( corrs_geo.positions );
  free( corrs_geo.colors_a );

}

void rendering::
init_hierarchical_corrs( bsc::gpu_geometry * geo,
                      const SequenceData * seq_data,
                      const ViewOptions * view_opts )
{
  // TODO: Add different coloring modes.
  if ( geo->vao != -1 ) 
  {
    bsc::free_gpu_geo( geo );
  }

  const std::unordered_map<FETShape*, Proxy*> * shape_to_proxy = 
                                                seq_data->model->shape_to_proxy;
  const std::vector<FETCorrespondence*> & corrs = seq_data->hierarchical_corrs;
  i32 n_corrs = corrs.size(); 

  bsc::geometry_data corrs_geo;
  corrs_geo.positions  = (r32*)malloc( 6 * n_corrs * sizeof(r32)  );
  corrs_geo.colors_a   = (u8*)malloc( 8 * n_corrs * sizeof(u8) ); 
  corrs_geo.user_data  = (r32*)malloc( 2 * n_corrs * sizeof(r32) );
  corrs_geo.n_vertices = 0;
  i32 pos_idx = 0, col_idx = 0, usr_idx = 0;

  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        corr_idx++ )
  {
    FETCorrespondence * corr = corrs[ corr_idx ];
    if(!view_opts->show_relationships_by_type[corr->RelationshipType()]) continue;
    if (!corr) { printf("Struct corr invalid\n"); continue; }   
    FETFeature * f1 = corr->Feature(0);
    FETFeature * f2 = corr->Feature(1);
    if ( !f1 || !f2 ) { printf("Struct corr features invalid\n"); continue; }    
    FETShape * s1 = f1->Shape();
    FETShape * s2 = f2->Shape();
    if ( !s1 || !s2 ) { printf("Struct corr shapes invalid\n"); continue; }

    R3Point p1 = f1->Position();
    p1.Transform( s1->Transformation( view_opts->xform_scheme ) );
    R3Point p2 = f2->Position();
    p2.Transform( s2->Transformation( view_opts->xform_scheme ) );

    R3Vector n1 = f1->Normal();
    n1.Transform( s1->Transformation( view_opts->xform_scheme ) );
    n1 = -n1;
    R3Vector n2 = f2->Normal();
    n2.Transform( s2->Transformation( view_opts->xform_scheme ) );
    n2 = -n2;

    i32 l1 = 0;
    i32 l2 = 0;
    
    // TODO: Put in a separate function
    auto shape_proxy1 = shape_to_proxy->find( f1->Shape() );
    auto shape_proxy2 = shape_to_proxy->find( f2->Shape() );

    if ( shape_proxy1 != shape_to_proxy->end() ) 
    {
      l1 = shape_proxy1->second->level + 1;
    }

    if ( shape_proxy2 != shape_to_proxy->end() ) 
    {
      l2 = shape_proxy2->second->level + 1;
    }

    // if level is 0, then just check shape.
    // else get inlier shape ids, and check if they are contained.
    bool should_show = 0;
    if ( l1 == 0 && is_within( f1->Shape()->reconstruction_index, view_opts ) )
    {
      should_show |= 1;
    }
    if ( l2 == 0 && is_within( f2->Shape()->reconstruction_index, view_opts ) )
    {
      should_show |= 1;
    }
    if ( l1 > 0 && is_within( shape_proxy1->second, view_opts ) )
    {
      should_show |= 1;
    }
    if ( l2 > 0 && is_within( shape_proxy2->second, view_opts ) )
    {
      should_show |= 1;
    }
    if ( !should_show ) continue;


    if( view_opts->inspection_mode && 
        view_opts->selected_class_id == PROXY_ID &&
        shape_proxy1->second->structure_index != view_opts->selected_object_id &&
        shape_proxy2->second->structure_index != view_opts->selected_object_id )
    {
      continue;
    }


    p1 += l1 * view_opts->proxy_offset * n1;
    p2 += l2 * view_opts->proxy_offset * n2;

    bsc::vec3i col_A = rel_type_map[ corr->RelationshipType() % 4 ];

    // copy in data
    corrs_geo.positions[pos_idx++] = p1.X(); 
    corrs_geo.positions[pos_idx++] = p1.Y(); 
    corrs_geo.positions[pos_idx++] = p1.Z(); 
                                          
    corrs_geo.positions[pos_idx++] = p2.X(); 
    corrs_geo.positions[pos_idx++] = p2.Y(); 
    corrs_geo.positions[pos_idx++] = p2.Z(); 
                                          
    corrs_geo.colors_a[col_idx++]  = corr->affinity * col_A.r;
    corrs_geo.colors_a[col_idx++]  = corr->affinity * col_A.g;
    corrs_geo.colors_a[col_idx++]  = corr->affinity * col_A.b;
    corrs_geo.colors_a[col_idx++]  = 255;

    corrs_geo.colors_a[col_idx++]  = corr->affinity * col_A.r;
    corrs_geo.colors_a[col_idx++]  = corr->affinity * col_A.g;
    corrs_geo.colors_a[col_idx++]  = corr->affinity * col_A.b;
    corrs_geo.colors_a[col_idx++]  = 255;

    corrs_geo.user_data[usr_idx++] = corr_idx;
    corrs_geo.user_data[usr_idx++] = corr_idx;

    corrs_geo.n_vertices += 2;
  }

  if ( corrs_geo.n_vertices > 0 )
  {
    bsc::init_gpu_geo(&corrs_geo, geo, POSITION | COLOR_A | USER_DATA );
  }

  free( corrs_geo.positions );
  free( corrs_geo.colors_a );
  free( corrs_geo.user_data );

}

void rendering::
init_cp_corrs( bsc::gpu_geometry * geo,
               const SequenceData * seq_data,
               const ViewOptions * view_opts )
{
  if ( geo->vao != -1 ) 
  {
    bsc::free_gpu_geo( geo );
  }

  const std::vector<FETCorrespondence*> & corrs = seq_data->cp_corrs;
  i32 n_corrs = corrs.size(); 

  bsc::geometry_data corrs_geo;
  corrs_geo.positions  = (r32*)malloc( 6 * n_corrs * sizeof(r32)  );
  corrs_geo.colors_a   = (u8*)malloc( 8 * n_corrs * sizeof(u8) ); 
  corrs_geo.n_vertices = 0;
  i32 pos_idx = 0, col_idx = 0;

  i32 lim_low  = view_opts->frame_idx_A;
  i32 lim_high = view_opts->frame_idx_B;

  for ( i32 corr_idx = 0 ;
        corr_idx < n_corrs ;
        corr_idx++ )
  {
    FETCorrespondence * corr = corrs[ corr_idx ];
    if (!corr) { continue; }   
    FETFeature * f1 = corr->Feature(0);
    FETFeature * f2 = corr->Feature(1);
    if (!f1 || !f2) { printf("CP corr features invalid\n"); continue; }   
    FETShape * s1 = f1->Shape();
    FETShape * s2 = f2->Shape();
    if (!s1 || !s2) { printf("CP corr features invalid\n"); continue; } 

    i32 idx_A = s1->reconstruction_index;
    i32 idx_B = s2->reconstruction_index;

    bool outside_range = idx_A < lim_low || idx_A > lim_high ||
                         idx_B < lim_low || idx_B > lim_high;
    bool wrong_idx = !( (idx_A == lim_low && idx_B == lim_high) ||
                        (idx_A == lim_high && idx_B == lim_low) );

    bool should_skip = view_opts->show_only_pair ? wrong_idx : outside_range;

    if ( should_skip ) continue;

    R3Point p1 = f1->Position();
    R3Point p2 = f2->Position();
    bsc::vec4d pos_A  = bsc::vec4d( p1.X(), p1.Y(), p1.Z(), 1.0f );
    bsc::vec4d pos_B  = bsc::vec4d( p2.X(), p2.Y(), p2.Z(), 1.0f );

    pos_A = view_opts->xforms_ptr[ idx_A ] * pos_A; 
    pos_B = view_opts->xforms_ptr[ idx_B ] * pos_B; 

    bsc::vec3i col_A = qual_map[ f1->GeneratorType() % 9 ];
    bsc::vec3i col_B = qual_map[ f1->GeneratorType() % 9 ];

    // copy in data
    corrs_geo.positions[pos_idx++] = pos_A.x; 
    corrs_geo.positions[pos_idx++] = pos_A.y; 
    corrs_geo.positions[pos_idx++] = pos_A.z; 
                                           
    corrs_geo.positions[pos_idx++] = pos_B.x; 
    corrs_geo.positions[pos_idx++] = pos_B.y; 
    corrs_geo.positions[pos_idx++] = pos_B.z; 
                                           
    corrs_geo.colors_a[col_idx++]    = corr->affinity*col_A.r;
    corrs_geo.colors_a[col_idx++]    = corr->affinity*col_A.g;
    corrs_geo.colors_a[col_idx++]    = corr->affinity*col_A.b;
    corrs_geo.colors_a[col_idx++]    = 255;

    corrs_geo.colors_a[col_idx++]    = corr->affinity*col_B.r;
    corrs_geo.colors_a[col_idx++]    = corr->affinity*col_B.g;
    corrs_geo.colors_a[col_idx++]    = corr->affinity*col_B.b;
    corrs_geo.colors_a[col_idx++]    = 255;

    corrs_geo.n_vertices += 2;
  }

  if ( corrs_geo.n_vertices > 0 )
  {
    bsc::init_gpu_geo(&corrs_geo, geo, POSITION | COLOR_A);
  }

  free( corrs_geo.positions );
  free( corrs_geo.colors_a );

}

void rendering::
init_frustum( bsc::gpu_geometry *cmd, const bsc::mat3 *K )
{
  i32 w = (*K)[2][0] * 2.0f;
  i32 h = (*K)[2][1] * 2.0f;

  bsc::vec2 coords[4] = {bsc::vec2(0, 0),
                         bsc::vec2(0, h),
                         bsc::vec2(w, h),
                         bsc::vec2(w, 0)};
  r32 depths[2] = {0.01, 0.2};
  bsc::vec3 tmp_positions[8];

  for (int j = 0; j < 2; ++j)
  {
    r32 depth = depths[j];
    for (int i = 0; i < 4; ++i)
    {
      i32 cx = coords[i].x;
      i32 cy = coords[i].y;

      bsc::vec3 pos;
      pos.x = ((cx + 0.5) - (*K)[2][0]) * depth / (*K)[0][0];
      pos.y = ((cy + 0.5) - (*K)[2][1]) * depth / (*K)[1][1];
      pos.z = -depth;

      tmp_positions[j * 4 + i] = pos;
    }
  }

  bsc::vec3 positions[16] = {tmp_positions[0], tmp_positions[1],
                             tmp_positions[2], tmp_positions[3],
                             tmp_positions[0], tmp_positions[4],
                             tmp_positions[5], tmp_positions[1],
                             tmp_positions[5], tmp_positions[6],
                             tmp_positions[2], tmp_positions[6],
                             tmp_positions[7], tmp_positions[3],
                             tmp_positions[7], tmp_positions[4]};

  bsc::geometry_data frustum;
  frustum.positions = (float *)&(positions[0]);
  frustum.n_vertices = 16;

  bsc::init_gpu_geo(&frustum, cmd, POSITION);
}

void rendering::
init_grid(bsc::gpu_geometry *cmd, i32 grid_size, r32 grid_spacing )
{
  using namespace bsc;
  vec3 *grid_locs = (vec3 *)malloc((grid_size + 1) * 4 * sizeof(vec3));

  r32 max_length = (grid_size * grid_spacing / 2.0f);
  i32 loc_idx = 0;
  for (i32 i = -grid_size / 2; i <= grid_size / 2; i++)
  {
    grid_locs[loc_idx++] = vec3(i * grid_spacing, 0.0f, -max_length);
    grid_locs[loc_idx++] = vec3(i * grid_spacing, 0.0f, max_length);
    grid_locs[loc_idx++] = vec3(-max_length, 0.0f, i * grid_spacing);
    grid_locs[loc_idx++] = vec3(max_length, 0.0f, i * grid_spacing);
  }

  geometry_data grid;
  grid.positions = (r32 *)&(grid_locs[0]);
  grid.n_vertices = (grid_size + 1) * 4;
  init_gpu_geo( &grid, cmd, POSITION);

  free(grid_locs);
}

void rendering::
init_disc(bsc::gpu_geometry *gpu_geo, i32 divs, r32 radius )
{
  using namespace bsc;

  vec3 *disc_pos = (vec3 *)malloc( (divs + 2) * sizeof(vec3) );

  i32 loc_idx = 0;
  disc_pos[loc_idx++] = vec3f( 0.0f, 0.0f, 0.0f );

  for ( r32 i = 0.0f; i <= 2.0f * RN_PI ; i += 2.0f * RN_PI / divs )
  {
    vec3 pos( radius * cos( i ), 0.0f, radius * sin( i ) );
    disc_pos[loc_idx++] = pos;
  }
  disc_pos[loc_idx++] = vec3f( radius * cos( 0.0f ), 0.0f, radius * sin( 0.0f ) );

  geometry_data disk_geo;
  disk_geo.positions = (r32 *)&(disc_pos[0]);
  disk_geo.n_vertices = divs + 2;
  init_gpu_geo( &disk_geo, gpu_geo, POSITION);

  free(disc_pos);
}

void rendering::
frame_to_pointcloud( bsc::geometry_data * pointcloud,
                     const bsc::img *color_img,
                     const bsc::img_u16 *depth_img,
                     const bsc::mat3 K,
                     const i32 downsample_factor,
                     const bool bitshift,
                     const r32 near_dist,
                     const r32 far_dist )
{

  // TODO(maciej): Normal estimation + filtering!!

  i32 valid_i = 0;
  i32 w = depth_img->width;
  i32 h = depth_img->height;

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
      if ( x == 0 || x == w - 1 ||
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
  // grd.BilateralFilter(3.0, 0.05);
 
  w = grd.XResolution();
  h = grd.YResolution();
  bsc::mat3 L( K );
  L[2][0] = (1.0f / downsample_factor) * K[2][0];
  L[2][1] = (1.0f / downsample_factor) * K[2][1];
  L[0][0] = (1.0f / downsample_factor) * K[0][0];
  L[1][1] = (1.0f / downsample_factor) * K[1][1];

  bsc::img_u16 depth_final( w, h, 1 );
  for ( int i = 0 ; i < w * h ; ++i )
  {
    (*depth_final.at(i)) = grd.GridValue(i) * 1000.0f;
  }

  r32 ratio_x = color_img->width / (r32)w; 
  r32 ratio_y = color_img->height / (r32)h; 
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

      bsc::vec3 p0( ((x) - L[2][0]) * depth0 / L[0][0],
                    ((h - y) - L[2][1]) * depth0 / L[1][1],
                    -depth0 );
      bsc::vec3 p1( ((x + 1) - L[2][0]) * depth1 / L[0][0],
                    ((h - y) - L[2][1]) * depth1 / L[1][1],
                     -depth1 );
      bsc::vec3 p2( ((x) - L[2][0]) * depth2 / L[0][0],
                    ((h - y - 1) - L[2][1]) * depth2 / L[1][1],
                    -depth2 );

      bsc::vec3 v0 = bsc::normalize( p1 - p0 );
      bsc::vec3 v1 = bsc::normalize( p2 - p0 );

      bsc::vec3 n = bsc::normalize( bsc::cross(v0, v1) );
      bsc::vec3 ray = bsc::normalize( bsc::vec3( 0.0f, 0.0f, 0.0f ) - p0 );
      if ( bsc::dot(n, ray) < 0 ) n = -n; 

      pointcloud->positions[3 * valid_i + 0] = p0.x;
      pointcloud->positions[3 * valid_i + 1] = p0.y;
      pointcloud->positions[3 * valid_i + 2] = p0.z;

      pointcloud->normals[3 * valid_i + 0] = n.x;
      pointcloud->normals[3 * valid_i + 1] = n.y;
      pointcloud->normals[3 * valid_i + 2] = n.z;

      pointcloud->colors_a[4 * valid_i + 0] = color[0];
      pointcloud->colors_a[4 * valid_i + 1] = color[1];
      pointcloud->colors_a[4 * valid_i + 2] = color[2];
      pointcloud->colors_a[4 * valid_i + 3] = 255;

      pointcloud->colors_b[4 * valid_i + 0] = color[0];
      pointcloud->colors_b[4 * valid_i + 1] = color[1];
      pointcloud->colors_b[4 * valid_i + 2] = color[2];
      pointcloud->colors_b[4 * valid_i + 3] = 255;

      valid_i++;
    }
  }

  pointcloud->n_vertices = valid_i;
  // exit(-1);
}


i32 rendering::
init_pointcloud( bsc::gpu_geometry *cmd,
                 const bsc::img *color_img,
                 const bsc::img_u16 *depth_img,
                 const bsc::mat3 K,
                 const r32 downsample_factor,
                 const bool bitshift,
                 const r32 near_dist,
                 const r32 far_dist )
{
  if ( cmd->vao <= 0 )
  {
    bsc::free_gpu_geo ( cmd );
  }
  
  bsc::geometry_data pointcloud;
  pointcloud.n_vertices = depth_img->width * depth_img->height;
  if ( pointcloud.n_vertices )
  {
    pointcloud.positions = (r32*)calloc(3 * pointcloud.n_vertices, sizeof(r32));
    pointcloud.normals   = (r32*)calloc(3 * pointcloud.n_vertices, sizeof(r32));
    pointcloud.colors_a  = (u8*)calloc(4 * pointcloud.n_vertices, sizeof(u8));
    pointcloud.colors_b  = (u8*)calloc(4 * pointcloud.n_vertices, sizeof(u8));
  
    frame_to_pointcloud( &pointcloud, color_img, depth_img, K, 
                        downsample_factor, bitshift, near_dist, far_dist);

    bsc::init_gpu_geo( &pointcloud, cmd, POSITION | NORMAL | COLOR_A | COLOR_B );

    free( pointcloud.positions );
    free( pointcloud.normals );
    free( pointcloud.colors_a );
    free( pointcloud.colors_b );
  }
  return pointcloud.n_vertices;
}


#endif //RENDERING_IMPLEMENTATION
#endif //_RENDERING_H_