// Source file for the rgbd loader program



////////////////////////////////////////////////////////////////////////
// Include files 
////////////////////////////////////////////////////////////////////////

#include "FET/FET.h"
#include "R3Surfels/R3Surfels.h"



////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

// File input/output options

static const char *input_configuration_name = NULL;
static const char *input_depth_directory = NULL;
static const char *input_color_directory = NULL;
static const char *output_reconstruction_name = NULL;
static const char *output_ssa_name = NULL;
static const char *output_ssb_name = NULL;
static const char *output_mesh_directory = NULL;


// Feature/surfel creation options

static RNBoolean create_boundary_features = FALSE;
static RNBoolean create_uniform_features = FALSE;
static RNBoolean create_plane_features = FALSE;
static RNBoolean create_crease_features = FALSE;
static RNBoolean create_multiresolution_surfels = FALSE;


// Image/pixel selection options

static int load_images_starting_at_index = 0;
static int load_images_ending_at_index = INT_MAX;
static int load_every_kth_image = 1;
static int max_image_resolution = 0;
static double min_feature_spacing = RN_EPSILON;
static double max_depth = 0;
static int omit_corners = 0;


// Algorithmic options

static double depth_blur_xy_sigma = 3;
static double depth_blur_d_sigma = 0.05;
static double normal_neighborhood_world_radius = 0.25;
static int normal_neighborhood_pixel_radius = 8; // could be 4 with less quality, could be 16 with better quality
static int normal_neighborhood_search = TRUE;
static int normal_ransac_iterations = 0;


// Printing options

static int print_verbose = 0;
static int print_debug = 0;



////////////////////////////////////////////////////////////////////////
// SUPPORT FILES
////////////////////////////////////////////////////////////////////////

#include "normal.cpp"
#include "segmentation.cpp"



////////////////////////////////////////////////////////////////////////
// SURFEL I/O STUFF
////////////////////////////////////////////////////////////////////////

static R3SurfelScene *
OpenSurfelScene(const char *scene_name, const char *database_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate surfel scene
  R3SurfelScene *scene = new R3SurfelScene();
  if (!scene) {
    fprintf(stderr, "Unable to allocate scene\n");
    return NULL;
  }

  // Open surfel scene files
  if (!scene->OpenFile(scene_name, database_name, "w", "w")) {
    delete scene;
    return NULL;
  }

  // Create surfel features 
  if (!CreateFeatures(scene)) exit(-1);

  // Print statistics
  if (print_verbose) {
    printf("Opened scene ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene->NObjects());
    printf("  # Labels = %d\n", scene->NLabels());
    printf("  # Assignments = %d\n", scene->NLabelAssignments());
    printf("  # Features = %d\n", scene->NFeatures());
    printf("  # Nodes = %d\n", scene->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene->Tree()->Database()->NBlocks());
    printf("  # Surfels = %d\n", scene->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Return scene
  return scene;
}



static int
CloseSurfelScene(R3SurfelScene *scene)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Print statistics
  if (print_verbose) {
    printf("Closing scene ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Objects = %d\n", scene->NObjects());
    printf("  # Labels = %d\n", scene->NLabels());
    printf("  # Assignments = %d\n", scene->NLabelAssignments());
    printf("  # Features = %d\n", scene->NFeatures());
    printf("  # Nodes = %d\n", scene->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene->Tree()->Database()->NBlocks());
    printf("  # Surfels = %d\n", scene->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Close surfel scene files
  if (!scene->CloseFile()) {
    delete scene;
    return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// INPUT STUFF
////////////////////////////////////////////////////////////////////////

static int
ReadIntrinsicsMatrix(R3Matrix& matrix, const char *filename, const char *dataset_format = NULL)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open matrix file %s\n", filename);
    return 0;
  }

  // Check file format
  if (dataset_format && !strcmp(dataset_format, "matterport")) {
    // Read intrinsic parameters
    double width, height, fx, fy, cx, cy, k1, k2, p1, p2, k3;
    if (fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &width, &height,
      &fx, &fy, &cx, &cy, &k1, &k2, &p1, &p2, &k3) != (unsigned int) 11) {
      fprintf(stderr, "Unable to read Matterport intrinsics matrix.\n");
      return 0;
    }

    // Fill in intrinsics matrix
    matrix = R3Matrix(fx, 0, cx,   0, fy, cy,   0, 0, 1);
  }
  else {
    // Read 3x3 matrix
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        fscanf(fp, "%lf", &matrix[i][j]);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



static int
ReadDepthImage(R2Grid& grid, const char *filename, const char *dataset_format = NULL)
{
  // Read depth image
  if (!grid.ReadFile(filename)) return 0;

  // Shift 3 bits (to compensate for shift done by SUNRGBD capture)
  if (dataset_format && (!strcmp(dataset_format, "sunrgbd") || !strcmp(dataset_format, "sun3d") || !strcmp(dataset_format, "nyu"))) {
    for (int i = 0; i < grid.NEntries(); i++) {
      unsigned int d = (unsigned int) (grid.GridValue(i) + 0.5);
      d = ((d >> 3) & 0x1FFF) | ((d & 0x7) << 13);
      grid.SetGridValue(i, d);
    }
  }

  // Convert integer coordinates into depths in meters
  if (strstr(filename, ".png")) {
    if (dataset_format && !strcmp(dataset_format, "matterport")) grid.Multiply(0.00025);
    else if (dataset_format && !strcmp(dataset_format, "tum")) grid.Multiply(0.0002);
    else if (dataset_format && (!strcmp(dataset_format, "sunrgbd") || !strcmp(dataset_format, "sun3d"))) grid.Multiply(0.001);
    else if (dataset_format && (!strcmp(dataset_format, "scannet") || !strcmp(dataset_format, "nyu"))) grid.Multiply(0.001);
    else grid.Multiply(0.001);
  }

#if 0
  // Process depth image 
  if (dataset_format && !strcmp(dataset_format, "matterport")) {
    const int pixel_radius = 4;
    const int min_pixel_count = 6; 
    const double max_delta_d = 0.2;
    const double max_closest_d = 0.025;

    // Fill in missing depth values
    if ((max_delta_d > 0) && (max_closest_d > 0)) {
      R2Grid copy(grid);
      for (int ix = 0; ix < grid.XResolution(); ix++) {
        for (int iy = 0; iy < grid.YResolution(); iy++) {
          RNScalar depth = copy.GridValue(ix, iy);
          if (RNIsPositive(depth)) continue;

          // Gather statistics from neighborhood
          int count = 0;
          RNScalar sum = 0;
          int closest_dd = INT_MAX;
          RNScalar closest_d = 0;
          RNScalar min_d = FLT_MAX;
          RNScalar max_d = -FLT_MAX;
          for (int s = -pixel_radius; s <= pixel_radius; s++) {
            if ((ix + s < 0) || (ix + s >= grid.XResolution())) continue;
            for (int t = -pixel_radius; t <= pixel_radius; t++) {
              if ((iy + t < 0) || (iy + t >= grid.YResolution())) continue;
              RNScalar d = copy.GridValue(ix+s, iy+t);
              if (RNIsNegativeOrZero(d)) continue;
              int dd = s*s + t*t;
              if (dd < closest_dd) { closest_dd = dd; closest_d = d; }
              if (d < min_d) min_d = d;
              if (d > max_d) max_d = d;
              sum += d;
              count++;
            }
          }

          // Fill in missing depth value with average if on planar surface
          if (count >= min_pixel_count) {
            if ((max_d - min_d) < max_delta_d) {
              RNScalar mean = sum / count;
              if (RNIsEqual(closest_d, mean, max_closest_d)) {
                grid.SetGridValue(ix, iy, mean);
              }
            }
          }
        }
      }
    }
  }
#endif
  
  // Smooth depth image
  if (dataset_format && (!strcmp(dataset_format, "sunrgbd") || !strcmp(dataset_format, "sun3d") ||
    !strcmp(dataset_format, "nyu") || !strcmp(dataset_format, "princeton"))) {
    if ((depth_blur_xy_sigma > 0) || (depth_blur_d_sigma > 0)) {
      RNScalar xy_sigma = depth_blur_xy_sigma * grid.XResolution() / 640.0;
      grid.Substitute(0, R2_GRID_UNKNOWN_VALUE);
      grid.BilateralFilter(xy_sigma, depth_blur_d_sigma);
      grid.Substitute(R2_GRID_UNKNOWN_VALUE, 0);
    }
  }

  // Print message
  if (0) {
    printf("Read depth image %s ...\n", filename);
    printf("  Resolution = %d %d\n", grid.XResolution(), grid.YResolution());
    printf("  Cardinality = %d\n", grid.Cardinality());
    printf("  Minimum = %g\n", grid.Minimum());
    printf("  Maximum = %g\n", grid.Maximum());
    printf("  L1Norm = %g\n", grid.L1Norm());
    printf("  L2Norm = %g\n", grid.L2Norm());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadColorImage(R2Image& image, const char *filename, const char *dataset_format = NULL)
{
  // Read color image
  if (!image.Read(filename)) return 0;

  // Print message
  if (0) {
    printf("Read color image %s ...\n", filename);
    printf("  Resolution = %d %d\n", image.Width(), image.Height());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// OUTPUT STUFF
////////////////////////////////////////////////////////////////////////

static int 
WriteReconstruction(FETReconstruction *reconstruction, const char *filename) 
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write file
  if (!reconstruction->WriteFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote reconstruction to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Shapes = %d\n", reconstruction->NShapes());
    printf("  # Matches = %d\n", reconstruction->NMatches());
    printf("  # Features = %d\n", reconstruction->NFeatures());
    printf("  # Correspondences = %d\n", reconstruction->NCorrespondences());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
WriteMesh(const char *filename, Segmentation *segmentation, const R2Grid& depth_image, 
  const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& radius_image, const R2Grid& boundary_image, const R2Image& color_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up)
{
  // Useful variables
  static const RNFlags bad_combination = R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG | R3_SURFEL_SHADOW_BOUNDARY_FLAG;

  // Just checking
  if (depth_image.XResolution() == 0) return 1;
  if (depth_image.YResolution() == 0) return 1;

  // Allocate temporary data
  R3Mesh *mesh = new R3Mesh();
  int max_vertices = depth_image.XResolution() * depth_image.YResolution();
  R3MeshVertex *vertex_block = new R3MeshVertex [ max_vertices ];
  R3MeshFace *face_block = new R3MeshFace [ 2 * max_vertices ];

  // Create vertices
  for (int j = 0; j < depth_image.YResolution(); j++) {
    for (int i = 0; i < depth_image.XResolution(); i++) {
      // Get vertex index
      int vertex_index = j*depth_image.XResolution() + i;

      // Check depth
      RNScalar depth = depth_image.GridValue(i, j);
      if (RNIsNegativeOrZero(depth)) continue;
      if ((max_depth > 0) && (depth > max_depth)) continue;

      // Get position
      RNCoord px = px_image.GridValue(i, j);
      RNCoord py = py_image.GridValue(i, j);
      RNCoord pz = pz_image.GridValue(i, j);
      R3Point position(px, py, pz);

      // Get normal
      RNCoord nx = nx_image.GridValue(i, j);
      RNCoord ny = ny_image.GridValue(i, j);
      RNCoord nz = nz_image.GridValue(i, j);
      R3Vector normal(nx, ny, nz);

      // Get color
      RNRgb color(0.0, 0.8, 0.0);
      if ((i < color_image.Width()) && (j < color_image.Height())) color = color_image.PixelRGB(i, j);

      // Create vertex
      mesh->CreateVertex(position, normal, color, &vertex_block[vertex_index]);
    }
  }

  // Create faces
  for (int i = 0; i < depth_image.XResolution()-1; i++) {
    for (int j = 0; j < depth_image.YResolution()-1; j++) {
      int face_index = 2*(j*depth_image.XResolution() + i);

      // Get depths
      RNScalar d00 = depth_image.GridValue(i, j);
      RNScalar d10 = depth_image.GridValue(i+1, j);
      RNScalar d01 = depth_image.GridValue(i, j+1);
      RNScalar d11 = depth_image.GridValue(i+1, j+1);

      // Get boundary flags
      RNFlags b00 = (unsigned int) (boundary_image.GridValue(i, j) + 0.5);
      RNFlags b10 = (unsigned int) (boundary_image.GridValue(i+1, j) + 0.5);
      RNFlags b01 = (unsigned int) (boundary_image.GridValue(i, j+1) + 0.5);
      RNFlags b11 = (unsigned int) (boundary_image.GridValue(i+1, j+1) + 0.5);

      // Get vertices
      R3MeshVertex *v00 = &vertex_block[(j+0)*depth_image.XResolution() + (i+0)];
      R3MeshVertex *v01 = &vertex_block[(j+1)*depth_image.XResolution() + (i+0)];
      R3MeshVertex *v10 = &vertex_block[(j+0)*depth_image.XResolution() + (i+1)];
      R3MeshVertex *v11 = &vertex_block[(j+1)*depth_image.XResolution() + (i+1)];

      // Create face1
      if (RNIsPositive(d00) && RNIsPositive(d10) && RNIsPositive(d11)) {
        if ((b00 | b10 | b11) != bad_combination) {
          mesh->CreateFace(v00, v10, v11, &face_block[face_index+0]);
        }
      } 

      // Create face2
      if (RNIsPositive(d00) && RNIsPositive(d11) && RNIsPositive(d01)) {
        if ((b00 | b11 | b01) != bad_combination) {
          mesh->CreateFace(v00, v11, v01, &face_block[face_index+1]);
        }
      } 
    }
  }

  // Write mesh
  mesh->WriteFile(filename);

  // Delete temporary data
  delete mesh;
  delete [] face_block;
  delete [] vertex_block;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// IMAGE RESAMPLING STUFF
////////////////////////////////////////////////////////////////////////

static int
ResampleDepthImage(R2Grid& image, R3Matrix& intrinsics_matrix, int xres, int yres)
{
  // Just checking
  if ((xres == 0) || (yres == 0)) return 0;
  if ((image.XResolution() == 0) || (image.YResolution() == 0)) return 0;

  // Compute scale factors
  double xscale = (double) image.XResolution() / (double) xres;
  double yscale = (double) image.YResolution() / (double) yres;

  // Copy original image and then initialize to zeroes
  R2Grid copy_image(image);
  image = R2Grid(xres, yres);

  // Resample images
  for (int iy = 0; iy < yres; iy++) {
    double cy = yscale * iy;
    int min_y = cy - 0.5*yscale + 0.5;
    if (min_y >= copy_image.YResolution()) continue;
    int max_y = cy + 0.5*yscale + 0.5;
    if (max_y < 0) continue;
    if (min_y < 0) min_y = 0;
    if (max_y >= copy_image.YResolution()) max_y = copy_image.YResolution()-1;
    if (max_y < min_y) max_y = min_y;
    
    for (int ix = 0; ix < xres; ix++) {
      double cx = xscale * ix;
      int min_x = cx - 0.5*xscale + 0.5;
      if (min_x >= copy_image.XResolution()) continue;
      int max_x = cx + 0.5*xscale + 0.5;
      if (max_x < 0) continue;
      if (min_x < 0) min_x = 0;
      if (max_x >= copy_image.XResolution()) max_x = copy_image.XResolution()-1;
      if (max_x < min_x) max_x = min_x;

      // Find minimum depth in neighborhood
      RNScalar min_depth = FLT_MAX;
      for (int j = min_y; j <= max_y; j++) {
        for (int i = min_x; i <= max_x; i++) {
          RNScalar depth = copy_image.GridValue(i, j);
          if ((depth <= 0) || (depth == R2_GRID_UNKNOWN_VALUE)) continue;
          if (depth < min_depth) min_depth = depth;
        }
      }

      // Set depth
      if (min_depth < FLT_MAX) {
        image.SetGridValue(ix, iy, min_depth);
      }
    }
  }

  // Update depth intrinsics matrix
  intrinsics_matrix[0][0] /= xscale;
  intrinsics_matrix[0][2] /= xscale;
  intrinsics_matrix[1][1] /= yscale;
  intrinsics_matrix[1][2] /= yscale;

  // Return success
  return 1;
}



static int
ResampleColorImage(R2Image& image, R3Matrix& intrinsics_matrix, int xres, int yres)
{
  // Just checking
  if ((xres == 0) || (yres == 0)) return 0;
  if ((image.Width() == 0) || (image.Height() == 0)) return 0;

  // Compute scale factors
  double xscale = (double) image.Width() / (double) xres;
  double yscale = (double) image.Height() / (double) yres;
  if (RNIsNotEqual(xscale, yscale, 0.01)) {
    fprintf(stderr, "Warning: anisotropic scaling of color image by factor %g\n", yscale/xscale);
  }
  
  // Copy original image and then initialize to zeroes
  R2Image copy_image(image);
  image = R2Image(xres, yres, 3);

  // Resample images
  for (int iy = 0; iy < yres; iy++) {
    double cy = yscale * iy;
    int min_y = cy - 0.5*yscale + 0.5;
    if (min_y >= copy_image.Height()) continue;
    int max_y = cy + 0.5*yscale + 0.5;
    if (max_y < 0) continue;
    if (min_y < 0) min_y = 0;
    if (max_y >= copy_image.Height()) max_y = copy_image.Height()-1;
    if (max_y < min_y) max_y = min_y;

    for (int ix = 0; ix < xres; ix++) {
      double cx = xscale * ix;
      int min_x = cx - 0.5*xscale + 0.5;
      if (min_x >= copy_image.Width()) continue;
      int max_x = cx + 0.5*xscale + 0.5;
      if (max_x < 0) continue;
      if (min_x < 0) min_x = 0;
      if (max_x >= copy_image.Width()) max_x = copy_image.Width()-1;
      if (max_x < min_x) max_x = min_x;

      // Compute weighted sum of colors in neighborhood
      RNRgb sum(0,0,0);
      RNScalar weight = 0;
      for (int j = min_y; j <= max_y; j++) {
        for (int i = min_x; i <= max_x; i++) {
          RNRgb color = copy_image.PixelRGB(i, j);
          if (color == RNblack_rgb) continue; // ??? for borders
          RNScalar w = 1.0;
          sum += w * color;
          weight += w;
        }
      }

      // Compute average
      if (weight > 0) {
        RNRgb color = sum / weight;
        image.SetPixelRGB(ix, iy, color);
      }
    }
  }

  // Update depth intrinsics matrix
  intrinsics_matrix[0][0] /= xscale;
  intrinsics_matrix[0][2] /= xscale;
  intrinsics_matrix[1][1] /= yscale;
  intrinsics_matrix[1][2] /= yscale;

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// RGBD PROCESSING STUFF
////////////////////////////////////////////////////////////////////////

static int
CreatePositionImages(const R2Grid& depth_image, 
  const R3Matrix& depth_intrinsics_matrix, const R4Matrix& extrinsics_matrix,
  R2Grid& px, R2Grid& py, R2Grid& pz)
{
  // Initialize position images
  px = depth_image; 
  px.Clear(0);
  py = px;
  pz = px;

  // Fill position images
  for (int i = 0; i < depth_image.XResolution(); i++) {
    for (int j = 0; j < depth_image.YResolution(); j++) {
      // Get depth
      RNScalar depth = depth_image.GridValue(i, j);
      if (RNIsNegativeOrZero(depth)) continue;

      // Get position in camera coordinate system 
      RNScalar x = ((i+0.5) - depth_intrinsics_matrix[0][2]) * depth / depth_intrinsics_matrix[0][0];
      RNScalar y = ((j+0.5) - depth_intrinsics_matrix[1][2]) * depth / depth_intrinsics_matrix[1][1];
      RNScalar z = -depth;
      R3Point position(x, y, z);

      // Transform by extrinsics matrix
      position = extrinsics_matrix * position;

      // Fill position images
      px.SetGridValue(i, j, position.X());
      py.SetGridValue(i, j, position.Y());
      pz.SetGridValue(i, j, position.Z());
    }
  }

  // Return success 
  return 1;
}



static int
CreateBoundaryImage(const R2Grid& depth_image, R2Grid& boundary_image, RNScalar depth_threshold = 0.1)
{
  // Initialize position images
  boundary_image = depth_image; 
  boundary_image.Clear(0);

  // Mark bottom and top border boundaries
  for (int i = 0; i < boundary_image.XResolution(); i++) {
    for (int j = 0; j < boundary_image.YResolution(); j++) {
      boundary_image.SetGridValue(i, j, R3_SURFEL_BORDER_BOUNDARY_FLAG);
      if (RNIsPositive(depth_image.GridValue(i, j))) break;
    }
    for (int j = boundary_image.YResolution()-1; j >= 0; j--) {
      boundary_image.SetGridValue(i, j, R3_SURFEL_BORDER_BOUNDARY_FLAG);
      if (RNIsPositive(depth_image.GridValue(i, j))) break;
    }
  }

  // Mark left and right border boundaries
  for (int j = 0; j < boundary_image.YResolution(); j++) {
    for (int i = 0; i < boundary_image.XResolution(); i++) {
      boundary_image.SetGridValue(i, j, R3_SURFEL_BORDER_BOUNDARY_FLAG);
      if (RNIsPositive(depth_image.GridValue(i, j))) break;
    }
    for (int i = boundary_image.XResolution()-1; i >= 0; i--) {
      boundary_image.SetGridValue(i, j, R3_SURFEL_BORDER_BOUNDARY_FLAG);
      if (RNIsPositive(depth_image.GridValue(i, j))) break;
    }
  }

  // Create copy of depth image with holes filled
  R2Grid filled_depth_image(depth_image);
  filled_depth_image.Substitute(0, R2_GRID_UNKNOWN_VALUE);
  filled_depth_image.FillHoles();

  // Mark interior boundaries, silhouettes, and shadows
  for (int i = 1; i < boundary_image.XResolution()-1; i++) {
    for (int j = 1; j < boundary_image.YResolution()-1; j++) {
      // Get original depth
      RNScalar depth = depth_image.GridValue(i, j);

      // Check if in hole
      if (RNIsNegativeOrZero(depth)) {
        boundary_image.SetGridValue(i, j, R3_SURFEL_BORDER_BOUNDARY_FLAG);
        continue;
      }

      // Get filled depth
      depth = filled_depth_image.GridValue(i, j);
      if (RNIsNegativeOrZero(depth)) continue;

      // Check depth relative to horizontal neighbors
      for (int k = 0; k < 4; k++) {
        int s = (k < 3) ? -1 : 0;
        int t = (k < 3) ? k-1 : -1;

        // Get depth on one side
        RNScalar depthA = filled_depth_image.GridValue(i-s, j-t);
        if (RNIsNegativeOrZero(depthA)) {
          boundary_image.SetGridValue(i, j, R3_SURFEL_BORDER_BOUNDARY_FLAG);
          break;
        }

        // Get depth on other side
        RNScalar depthB = filled_depth_image.GridValue(i+s, j+t);
        if (RNIsNegativeOrZero(depthB)) {
          boundary_image.SetGridValue(i, j, R3_SURFEL_BORDER_BOUNDARY_FLAG);
          break;
        }

        // Check differences of depth for shadow/silhouette
        RNScalar deltaA = depth - depthA;
        RNScalar deltaB = depthB - depth;
        RNScalar threshold = depth * depth_threshold;
        if (threshold < 0.1) threshold = 0.1;
        if (deltaA < -threshold) {
          if (deltaA < 4*deltaB) {
            boundary_image.SetGridValue(i-s, j-t, R3_SURFEL_SHADOW_BOUNDARY_FLAG);
            boundary_image.SetGridValue(i, j, R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG);
          }
        }
        else if (deltaA > threshold) {
          if (deltaA > 4*deltaB) {
            boundary_image.SetGridValue(i-s, j-t, R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG);
            boundary_image.SetGridValue(i, j, R3_SURFEL_SHADOW_BOUNDARY_FLAG);
          }
        }
        if (deltaB < -threshold) {
          if (deltaB < 4*deltaA) {
            boundary_image.SetGridValue(i+s, j+t, R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG);
            boundary_image.SetGridValue(i, j, R3_SURFEL_SHADOW_BOUNDARY_FLAG);
          }
        }
        else if (deltaB > threshold) {
          if (deltaB > 4*deltaA) {
            boundary_image.SetGridValue(i+s, j+t, R3_SURFEL_SHADOW_BOUNDARY_FLAG);
            boundary_image.SetGridValue(i, j, R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG);
          }
        }
      }
    }
  }

  // Return success 
  return 1;
}



static int
CreateNormalImages(const R2Grid& depth_image, 
  const R2Grid& px, const R2Grid& py, const R2Grid& pz, const R2Grid& boundary_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  R2Grid& nx, R2Grid& ny, R2Grid& nz, R2Grid& radius)
{
  // Initialize normal images
  nx = depth_image; 
  nx.Clear(0);
  ny = nx;
  nz = nx;
  radius = nx;

  // Allocate array of neighborhood points
  int neighborhood_pixel_radius = normal_neighborhood_pixel_radius * depth_image.XResolution() / 640.0 + 0.5;
  if (neighborhood_pixel_radius == 0) neighborhood_pixel_radius = 1;
  int neighborhood_pixel_radius_squared = neighborhood_pixel_radius * neighborhood_pixel_radius;
  RNScalar neighborhood_world_radius_squared = normal_neighborhood_world_radius * normal_neighborhood_world_radius;
  R3Point *neighborhood_points = new R3Point [ depth_image.NEntries() ];
  // R3Point *neighborhood_points = new R3Point [ neighborhood_diameter * neighborhood_diameter ];

  // Allocate search grid
  R2Grid search_grid(depth_image.XResolution(), depth_image.YResolution());
  const RNScalar *search_valuesp = search_grid.GridValues();
  search_grid.Clear(-1);

  // Fill normal images
  for (int i = 0; i < depth_image.XResolution(); i++) {
    for (int j = 0; j < depth_image.YResolution(); j++) {
      int neighborhood_npoints = 0;

      // Check depth
      RNScalar depth = depth_image.GridValue(i, j);
      if (RNIsNegativeOrZero(depth)) continue;
      if ((max_depth > 0) && (depth > max_depth)) continue;

      // Get image position
      R2Point image_position(i, j);

      // Get world position
      RNScalar x = px.GridValue(i, j);
      RNScalar y = py.GridValue(i, j);
      RNScalar z = pz.GridValue(i, j);
      R3Point world_position(x, y, z);

      // Create array of neighborhood points
      if (normal_neighborhood_search) {
        // Initialize stack
        int seed_index;
        RNArray<const RNScalar *> stack;
        search_grid.IndicesToIndex(i, j, seed_index);
        const RNScalar *seed_valuep = &search_valuesp[seed_index];
        search_grid.SetGridValue(seed_index, seed_index);
        stack.Insert(seed_valuep);

        // Depth first search, not extending beyond boundaries
        while (!stack.IsEmpty()) {
          // Pop current point off stack
          int ix, iy;
          const RNScalar *current_valuep = stack.Tail(); stack.RemoveTail();
          int current_index = current_valuep - search_valuesp;
          search_grid.IndexToIndices(current_index, ix, iy);

          // Get current info
          RNScalar current_x = px.GridValue(ix, iy);
          RNScalar current_y = py.GridValue(ix, iy);
          RNScalar current_z = pz.GridValue(ix, iy);
          R3Point current_position(current_x, current_y, current_z);
          RNScalar current_boundary_value = boundary_image.GridValue(ix, iy);
          int current_boundary = (int) (current_boundary_value + 0.5);
        
          // Add point to array
          neighborhood_points[neighborhood_npoints++] = current_position;

          // Add adjacent points to stack
          for (int nx = ix-1; nx <= ix+1; nx++) {
            if ((nx < 0) || (nx >= search_grid.XResolution())) continue;
            for (int ny = iy-1; ny <= iy+1; ny++) {
              if ((ny < 0) || (ny >= search_grid.YResolution())) continue;

              // Check if neighbor has already been visited
              RNScalar neighbor_value = search_grid.GridValue(nx, ny);
              if (RNIsEqual(neighbor_value, seed_index)) continue;
              search_grid.SetGridValue(nx, ny, seed_index);

              // Check neighbor depth
              RNScalar neighbor_depth = depth_image.GridValue(nx, ny);
              if (RNIsNegativeOrZero(neighbor_depth)) continue;
              if ((max_depth > 0) && (depth > max_depth)) continue;

              // Check neighbor image distance
              R2Point neighbor_image_position(nx, ny);
              if (R2SquaredDistance(image_position, neighbor_image_position) > neighborhood_pixel_radius_squared) continue;
            
              // Check neighbor world distance
              RNScalar neighbor_x = px.GridValue(nx, ny);
              RNScalar neighbor_y = py.GridValue(nx, ny);
              RNScalar neighbor_z = pz.GridValue(nx, ny);
              R3Point neighbor_world_position(neighbor_x, neighbor_y, neighbor_z);
              if (neighborhood_world_radius_squared > 0) {
                if (R3SquaredDistance(world_position, neighbor_world_position) > neighborhood_world_radius_squared) continue;
              }

              // Check if neighbor is across boundary
              int neighbor_boundary = (int) (boundary_image.GridValue(nx, ny) + 0.5);
              if ((current_boundary == R3_SURFEL_SHADOW_BOUNDARY_FLAG) && (neighbor_boundary == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG)) continue;
              if ((current_boundary == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG) && (neighbor_boundary == R3_SURFEL_SHADOW_BOUNDARY_FLAG)) continue;
        
              // Add neighbor to search
              int neighbor_index;
              search_grid.IndicesToIndex(nx, ny, neighbor_index);
              const RNScalar *neighbor_valuep = &search_valuesp[neighbor_index];
              stack.Insert(neighbor_valuep);
            }
          }
        }
      }
      else {
        // Add neighbor points to array
        for (int nx = i-neighborhood_pixel_radius; nx <= i+neighborhood_pixel_radius; nx++) {
          if ((nx < 0) || (nx >= depth_image.XResolution())) continue;
          for (int ny = j-neighborhood_pixel_radius; ny <= j+neighborhood_pixel_radius; ny++) {
            if ((ny < 0) || (ny >= depth_image.YResolution())) continue;
            
            // Check neighbor depth
            RNScalar neighbor_depth = depth_image.GridValue(nx, ny);
            if (RNIsNegativeOrZero(neighbor_depth)) continue;
            if ((max_depth > 0) && (depth > max_depth)) continue;

            // Check neighbor world position
            RNScalar neighbor_x = px.GridValue(nx, ny);
            RNScalar neighbor_y = py.GridValue(nx, ny);
            RNScalar neighbor_z = pz.GridValue(nx, ny);
            R3Point neighbor_world_position(neighbor_x, neighbor_y, neighbor_z);
            if (neighborhood_world_radius_squared > 0) {
              if (R3SquaredDistance(world_position, neighbor_world_position) > neighborhood_world_radius_squared) continue;
            }

            // Add neighbor point to array
            neighborhood_points[neighborhood_npoints++] = neighbor_world_position;
          }
        }
      }

      // Check number of neighbor points
      if (neighborhood_npoints < 3) continue;

      // Solve for normal
      RNScalar variances[3];
      R3Point centroid = R3Centroid(neighborhood_npoints, neighborhood_points);
      R3Triad triad = R3PrincipleAxes(centroid, neighborhood_npoints, neighborhood_points, NULL, variances);
      R3Vector normal = triad[2];

      // Flip normal to point towards camera
      R3Plane plane(world_position, normal);
      if (R3SignedDistance(plane, viewpoint) < 0) normal.Flip(); 

      // Refine normal with RANSAC
      if (normal_ransac_iterations > 0) {      
        RNScalar max_inlier_fraction = 0;
        RNScalar avg_inlier_fraction = 0;
        RefineNormalWithRansac(normal, world_position, neighborhood_points, neighborhood_npoints, 0.01, 16,
          &max_inlier_fraction, &avg_inlier_fraction);
        plane.Reset(world_position, normal);
        if (R3SignedDistance(plane, viewpoint) < 0) normal.Flip();
      }

      // Fill normal images
      nx.SetGridValue(i, j, normal.X());
      ny.SetGridValue(i, j, normal.Y());
      nz.SetGridValue(i, j, normal.Z());

      // Fill radius image
      RNScalar r = sqrt(variances[0]);
      if (neighborhood_pixel_radius > 1) r /= normal_neighborhood_pixel_radius;
      radius.SetGridValue(i, j, r);
    }
  }

  // Delete array of neighborhood points
  delete [] neighborhood_points;

  // Return success 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Surfel Processing
////////////////////////////////////////////////////////////////////////

static R3SurfelBlock *
LoadSurfels(R3SurfelScene *scene, const RNArray<Point *>& points,
  R3SurfelObject *parent_object, R3SurfelNode *parent_node,
  const char *node_name, const R2Grid& depth_image, const R3Point& origin)
{
  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  R3SurfelDatabase *database = tree->Database();
  int cx = depth_image.XResolution()/2;
  int cy = depth_image.YResolution()/2;
  double rx = cx;
  double ry = cy;
  
  // Allocate array of surfels
  int nsurfels = 0;
  R3Surfel *surfels = new R3Surfel [ points.NEntries() ];
  if (!surfels) {
    fprintf(stderr, "Unable to allocate surfels for %s\n", node_name);
    return NULL;
  }

  // Load surfels into array
  for (int i = 0; i < points.NEntries(); i++) {
    Point *point = points.Kth(i);
        
    // Check depth
    if (RNIsNegativeOrZero(point->depth)) continue;
    if ((max_depth > 0) && (point->depth > max_depth)) continue;

    // Get/check grid index
    int ix, iy;
    if (point->grid_index < 0) continue;
    depth_image.IndexToIndices(point->grid_index, ix, iy);

    // Check if outside elipse centered in image
    if (omit_corners) {
      int dx = (ix - cx) / rx;
      int dy = (iy - cy) / ry;
      if (dx*dx + dy*dy > 1) continue;
    }

    // Set position
    float px = (float) (point->position.X() - origin.X());
    float py = (float) (point->position.Y() - origin.Y());
    float pz = (float) (point->position.Z() - origin.Z());
    surfels[nsurfels].SetCoords(px, py, pz);

    // Set normal
    float nx = (float) point->normal.X();
    float ny = (float) point->normal.Y();
    float nz = (float) point->normal.Z();
    surfels[nsurfels].SetNormal(nx, ny, nz);

    // Set radius
    float radius = (float) point->radius;
    surfels[nsurfels].SetRadius(radius);

    // Set color
    surfels[nsurfels].SetColor(point->color);

    // Set boundary flags
    if (point->boundary == R3_SURFEL_BORDER_BOUNDARY_FLAG) surfels[nsurfels].SetBorderBoundary(TRUE);
    else if (point->boundary == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG) surfels[nsurfels].SetSilhouetteBoundary(TRUE);
    else if (point->boundary == R3_SURFEL_SHADOW_BOUNDARY_FLAG) surfels[nsurfels].SetShadowBoundary(TRUE);

    // Increment number of surfels
    nsurfels++;
  }

  // Create block from array of surfels
  R3SurfelBlock *block = new R3SurfelBlock(surfels, nsurfels, origin);
  if (!block) {
    fprintf(stderr, "Unable to allocate block\n");
    delete [] surfels;
    return NULL;
  }

  // Delete array of surfels
  delete [] surfels;

  // Insert block into database
  block->UpdateProperties();
  database->InsertBlock(block);

  // Create node
  R3SurfelNode *node = NULL;
  if (block && parent_node) {
    // Create node
    node = new R3SurfelNode(node_name);
    if (!node) {
      fprintf(stderr, "Unable to allocate node for %s\n", node_name);
      delete block;
      return NULL;
    }
            
    // Insert node into tree
    tree->InsertNode(node, parent_node);
    node->InsertBlock(block);
    node->UpdateProperties();
  }

  // Creat object
  R3SurfelObject *object = NULL;
  if (node && parent_object) {
    // Create object
    object = new R3SurfelObject(node_name);
    if (!object) {
      fprintf(stderr, "Unable to allocate object for %s\n", node_name);
      delete block;
      delete node;
      return NULL;
    }

    // Insert object
    scene->InsertObject(object, parent_object);
    object->InsertNode(node);
    object->UpdateProperties();

    // Create PCA object property
    R3SurfelObjectProperty *pca = new R3SurfelObjectProperty(R3_SURFEL_OBJECT_PCA_PROPERTY, object);
    scene->InsertObjectProperty(pca);
  }

  // Release block
  database->ReleaseBlock(block);

  // Return block
  return block;
}



static int
LoadSurfels(R3SurfelScene *scene, 
  const char *scan_name, Segmentation *segmentation, const R2Grid& depth_image, 
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up)
{
  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  R3SurfelObject *parent_object = scene->RootObject();
  R3SurfelNode *parent_node = tree->RootNode();
  int width = depth_image.XResolution();
  int height = depth_image.YResolution();
  
  // Create array of points outside any cluster
  RNArray<Point *> unclustered_points;
  for (int i = 0; i < segmentation->points.NEntries(); i++) {
    Point *point = segmentation->points.Kth(i);
    if (!point->cluster) unclustered_points.Insert(point);
  }

  // Create node
  char node_name[1024];
  sprintf(node_name, "SCAN:%s", scan_name);
  R3SurfelNode *scan_node = new R3SurfelNode(node_name);
  if (!scan_node) {
    fprintf(stderr, "Unable to allocate node for %s\n", scan_name);
    return 0;
  }
            
  // Insert node into tree
  tree->InsertNode(scan_node, parent_node);

  // Load surfels for points inside every cluster
  for (int i = 0; i < segmentation->clusters.NEntries(); i++) {
    Cluster *cluster = segmentation->clusters.Kth(i);
    if (cluster->points.NEntries() == 0) continue;
    
    // Create node name
    char node_name[1024];
    sprintf(node_name, "SCAN:%s:%d", scan_name, i);

    // Load surfels
    LoadSurfels(scene, cluster->points, parent_object, scan_node, node_name, depth_image, viewpoint);
  }

  // Load surfels for points outside any cluster
  if (unclustered_points.NEntries() > 0) {
    // Create node name
    char node_name[1024];
    sprintf(node_name, "SCAN:%s:unclustered", scan_name);

    // Load surfels
    LoadSurfels(scene, unclustered_points, NULL, scan_node, node_name, depth_image, viewpoint);
  }
            
  // Update node properties
  scan_node->UpdateProperties();

  // Create scan
  R3SurfelScan *scan = new R3SurfelScan(scan_name);
  if (!scan) {
    fprintf(stderr, "Unable to allocate scan\n");
    return 0;
  }

  // Assign scan properties
  scan->SetViewpoint(viewpoint);
  scan->SetOrientation(towards, up);
  scan->SetImageDimensions(width, height);
  scan->SetImageCenter(R2Point(0.5*width, 0.5*height));
  scan->SetNode(scan_node);
          
  // Insert scan
  scene->InsertScan(scan);

  // Return success
  return 1;
}



static int
LoadSurfels(R3SurfelScene *scene,
  const char *scan_name, Segmentation *segmentation, const R2Grid& depth_image, 
  const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& radius_image, const R2Grid& boundary_image, const R2Image& color_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up)
{
  // Get convenient variables
  R3SurfelTree *tree = scene->Tree();
  R3SurfelDatabase *database = tree->Database();
  R3SurfelNode *parent_node = tree->RootNode();
  int width = depth_image.XResolution();
  int height = depth_image.YResolution();
  int cx = depth_image.XResolution()/2;
  int cy = depth_image.YResolution()/2;
  double rx = cx;
  double ry = cy;
  
  // Allocate array of surfels
  int max_surfels = depth_image.NEntries();
  R3Surfel *surfels = new R3Surfel [ max_surfels ];
  if (!surfels) {
    fprintf(stderr, "Unable to allocate surfels for %s\n", scan_name);
    return 0;
  }

  // Load surfels into array
  int nsurfels = 0;
  for (int i = 0; i < depth_image.XResolution(); i++) {
    for (int j = 0; j < depth_image.YResolution(); j++) {
      // Check depth
      RNScalar depth = depth_image.GridValue(i, j);
      if (RNIsNegativeOrZero(depth)) continue;
      if ((max_depth > 0) && (depth > max_depth)) continue;

      // Check if outside elipse centered in image
      if (omit_corners) {
        int dx = (i - cx) / rx;
        int dy = (j - cy) / ry;
        if (dx*dx + dy*dy > 1) continue;
      }

      // Set position
      float px = (float) (px_image.GridValue(i, j) - viewpoint.X());
      float py = (float) (py_image.GridValue(i, j) - viewpoint.Y());
      float pz = (float) (pz_image.GridValue(i, j) - viewpoint.Z());
      surfels[nsurfels].SetCoords(px, py, pz);

      // Set normal
      float nx = (float) (nx_image.GridValue(i, j));
      float ny = (float) (ny_image.GridValue(i, j));
      float nz = (float) (nz_image.GridValue(i, j));
      surfels[nsurfels].SetNormal(nx, ny, nz);

      // Set radius
      float radius = (float) (radius_image.GridValue(i, j));
      surfels[nsurfels].SetRadius(radius);

      // Set color
      RNRgb color(0.0, 0.8, 0.0);
      if ((i < color_image.Width()) && (j < color_image.Height())) {
        color = color_image.PixelRGB(i, j);
        surfels[nsurfels].SetColor(color);
      }

      // Set boundary flags
      int boundary = (int) (boundary_image.GridValue(i, j) + 0.5);
      if (boundary == R3_SURFEL_BORDER_BOUNDARY_FLAG) surfels[nsurfels].SetBorderBoundary(TRUE);
      else if (boundary == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG) surfels[nsurfels].SetSilhouetteBoundary(TRUE);
      else if (boundary == R3_SURFEL_SHADOW_BOUNDARY_FLAG) surfels[nsurfels].SetShadowBoundary(TRUE);

      // Increment number of surfels
      nsurfels++;
    }
  }

  // Create block from array of surfels
  R3SurfelBlock *block = new R3SurfelBlock(surfels, nsurfels, viewpoint);
  if (!block) {
    fprintf(stderr, "Unable to allocate block\n");
    delete [] surfels;
    return 0;
  }

  // Delete array of surfels
  delete [] surfels;

  // Update block properties
  block->UpdateProperties();
  
  // Insert block into database
  database->InsertBlock(block);

  // Create node name
  char node_name[1024];
  sprintf(node_name, "SCAN:%s", scan_name);
  
  // Create node
  R3SurfelNode *node = new R3SurfelNode(node_name);
  if (!node) {
    fprintf(stderr, "Unable to allocate node\n");
    delete block;
    return 0;
  }
            
  // Insert node into tree
  tree->InsertNode(node, parent_node);

  // Insert block into node
  node->InsertBlock(block);
  
  // Update node properties
  node->UpdateProperties();

  // Release block
  database->ReleaseBlock(block);

  // Create scan
  R3SurfelScan *scan = new R3SurfelScan(node_name);
  if (!scan) {
    fprintf(stderr, "Unable to allocate scan\n");
    return 0;
  }

  // Assign scan properties
  scan->SetViewpoint(viewpoint);
  scan->SetOrientation(towards, up);
  scan->SetImageDimensions(width, height);
  scan->SetImageCenter(R2Point(0.5*width, 0.5*height));
  scan->SetNode(node);
          
  // Insert scan
  scene->InsertScan(scan);

  // Return success
  return 1;
}



static int
CreateMultiresolutionSurfels(R3SurfelScene *scene)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Creating hierarchy ...\n");
    fflush(stdout);
  }

  // Get surfel tree
  R3SurfelTree *tree = scene->Tree();
  if (!tree) {
    fprintf(stderr, "Surfel scene has no tree\n");
    return 0;
  }    

  // Split nodes
  int max_parts_per_node = 8;
  int max_blocks_per_node = 32;
  RNScalar max_node_complexity = 1024;
  RNScalar max_block_complexity = 1024;
  RNLength max_leaf_extent = 1.0;
  RNLength max_block_extent = 1.0;
  int max_levels = 64;
  tree->SplitNodes(max_parts_per_node, max_blocks_per_node, 
    max_node_complexity, max_block_complexity, 
    max_leaf_extent, max_block_extent, max_levels);

  // Create multiresolution blocks
  RNScalar multiresolution_factor = 0.25;
  tree->CreateMultiresolutionBlocks(multiresolution_factor, max_block_complexity);

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Nodes = %d\n", scene->Tree()->NNodes());
    printf("  # Blocks = %d\n", scene->Tree()->Database()->NBlocks());
    printf("  # Surfels = %d\n", scene->Tree()->Database()->NSurfels());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// FET CREATION STUFF
////////////////////////////////////////////////////////////////////////

static FETReconstruction *
CreateReconstruction(void)
{            
  // Create reconstruction
  FETReconstruction *reconstruction = new FETReconstruction();    
  if (!reconstruction) {
    fprintf(stderr, "Unable to allocate reconstruction\n");
    return NULL;
  }
    
  // Return reconstruction
  return reconstruction;
}



static FETShape *
CreateShape(FETReconstruction *reconstruction, const char *scan_name,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up)
{            
  // Get shape name
  char shape_name[4096];
  if (scan_name) sprintf(shape_name, "SCAN:%s", scan_name);
  else sprintf(shape_name, "SCAN:%d", reconstruction->NShapes());

  // Create shape
  FETShape *shape = new FETShape(reconstruction);    
  if (!shape) {
    fprintf(stderr, "Unable to allocate shape\n");
    return NULL;
  }
    
  // Set shape properties
  shape->SetName(shape_name);
  shape->SetOrigin(viewpoint);
  shape->SetViewpoint(viewpoint);
  shape->SetTowards(towards);
  shape->SetUp(up);

  // Return shape
  return shape;
}



static FETFeature *
CreateFeature(FETReconstruction *reconstruction, FETShape *shape,
  int ix, int iy, int shape_type, int generator_type, const FETDescriptor& descriptor,
  const R3Point& position, const R3Vector& direction, const R3Vector& normal, const RNRgb& color,
  RNScalar depth, RNScalar radius, RNScalar salience, unsigned int boundary, 
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  R3Kdtree<FETFeature *>& previous_features, RNBoolean update_previous_features = TRUE)
{
  // Check depth
  if (RNIsNegativeOrZero(depth)) return NULL;

  // Check salience
  static const RNScalar min_salience = 1E-6;
  static const RNScalar max_salience = 1.0;
  if ((min_salience > 0) && (salience < min_salience)) return NULL;
  if ((max_salience > 0) && (salience > max_salience)) salience = max_salience;

  // Check spacing
  if (min_feature_spacing > 0) {
    if (previous_features.FindAny(position, 0, min_feature_spacing)) return NULL;
  }
  
  // Check radius
  if (min_feature_spacing > 0) {
    RNLength min_radius = 0.5 * min_feature_spacing;
    if (radius < min_radius) radius = min_radius;
  }

  // Compute flags
  RNFlags flags = 0;
  if (boundary != 0) flags.Add(boundary);
  if (shape_type == POINT_FEATURE_SHAPE) flags.Add(FEATURE_IS_POINT);
  if (shape_type == LINE_FEATURE_SHAPE) flags.Add(FEATURE_IS_LINEAR);

  // Create feature
  FETFeature *feature = new FETFeature(reconstruction, shape_type, position, direction, normal, radius, descriptor, color, flags);
  if (!feature) {
    fprintf(stderr, "Unable to create feature\n");
    return NULL;
  }

  // Set salience
  feature->SetSalience(salience);

  // Set generator
  feature->generator_type = generator_type;
  
  // Insert feature
  shape->InsertFeature(feature);

  // Add to set for min_feature_spacing checks
  if ((min_feature_spacing > 0) && (update_previous_features)) {
    previous_features.InsertPoint(feature);
  }

  // Return feature
  return feature;
}



static FETFeature *
CreateFeature(FETReconstruction *reconstruction, FETShape *shape,
  int ix, int iy, int shape_type, int generator_type, const FETDescriptor& descriptor,
  const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image, 
  const R2Grid& boundary_image, const R2Image& color_image, RNScalar salience,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  R3Kdtree<FETFeature *>& previous_features, RNBoolean update_previous_features = TRUE)
{
  // Get/check depth
  RNScalar depth = depth_image.GridValue(ix, iy);
  if (RNIsNegativeOrZero(depth)) return NULL;

  // Get position
  RNScalar px = px_image.GridValue(ix, iy);
  RNScalar py = py_image.GridValue(ix, iy);
  RNScalar pz = pz_image.GridValue(ix, iy);
  R3Point position(px, py, pz);

  // Get normal
  RNScalar nx = nx_image.GridValue(ix, iy);
  RNScalar ny = ny_image.GridValue(ix, iy);
  RNScalar nz = nz_image.GridValue(ix, iy);
  R3Vector normal(nx, ny, nz);

  // Get radius
  RNScalar radius = radius_image.GridValue(ix, iy);

  // Get color
  RNRgb color(0.0, 0.8, 0.0);
  if ((ix < color_image.Width()) && (iy < color_image.Height())) {
    color = color_image.PixelRGB(ix, iy);
  }

  // Get boundary
  unsigned int boundary = 0;
  unsigned int boundary_flags = (unsigned int) (boundary_image.GridValue(ix, iy) + 0.5);
  if (boundary_flags == R3_SURFEL_BORDER_BOUNDARY_FLAG) boundary = FEATURE_IS_ON_BORDER_BOUNDARY;
  else if (boundary_flags == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG) boundary = FEATURE_IS_ON_SILHOUETTE_BOUNDARY;
  else if (boundary_flags == R3_SURFEL_SHADOW_BOUNDARY_FLAG) boundary = FEATURE_IS_ON_SHADOW_BOUNDARY;

  // Get direction
  R3Vector direction(0,0,0);
  if ((generator_type == SILHOUETTE_FEATURE_TYPE) || (generator_type == SHADOW_FEATURE_TYPE)) {
    // Get neighbors
    int nneighbors = 0;
    const int image_radius = 2;
    R3Point neighbors[(2*image_radius + 1)*(2*image_radius + 1)];
    for (int s = -image_radius; s <= image_radius; s++) {
      if ((ix+s < 0) || (ix+s >= px_image.XResolution())) continue;
      for (int t = -image_radius; t <= image_radius; t++) {
        if ((iy+t < 0) || (iy+t >= px_image.YResolution())) continue;
        int nb = (int) (boundary_image.GridValue(ix+s,iy+t) + 0.5);
        if ((generator_type == SILHOUETTE_FEATURE_TYPE) && (!(nb == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG))) continue;
        if ((generator_type == SHADOW_FEATURE_TYPE) && (!(nb == R3_SURFEL_SHADOW_BOUNDARY_FLAG))) continue;
        R3Point neighbor(px_image.GridValue(ix+s,iy+t), py_image.GridValue(ix+s,iy+t), pz_image.GridValue(ix+s,iy+t));
        neighbors[nneighbors++] = neighbor;
      }
    }

    // Compute direction from principle axis
    direction = normal % towards;
    if (nneighbors > 2) {
      R3Point centroid = R3Centroid(nneighbors, neighbors);
      R3Triad triad = R3PrincipleAxes(centroid, nneighbors, neighbors);
      direction = triad[0];
    }

    // Orient direction
    R3Vector right = towards % up;
    if (direction.Dot(right) < -0.1) direction.Flip();
    if (direction.Dot(up) < -0.1) direction.Flip();
    direction.Normalize();

    // Check direction
    if (RNIsZero(direction.Length())) return NULL;
  }

  // Create feature
  return CreateFeature(reconstruction, shape,
    ix, iy, shape_type, generator_type, descriptor,
    position, direction, normal, color,
    depth, radius, salience, boundary, 
    viewpoint, towards, up,
    previous_features, update_previous_features);
}



static int
CreateBoundaryFeatures(FETReconstruction *reconstruction, FETShape *shape,
  const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image, 
  const R2Grid& boundary_image, const R2Image& color_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  R3Kdtree<FETFeature *>& previous_features, RNBoolean update_previous_features = TRUE)
{
  // Create silhouette and shadow image
  R2Grid component_image(boundary_image.XResolution(), boundary_image.YResolution());
  for (int i = 0; i < boundary_image.NEntries(); i++) {
    int boundary = (int) (boundary_image.GridValue(i) + 0.5);
    if ((boundary == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG) || (boundary == R3_SURFEL_SHADOW_BOUNDARY_FLAG)) {
      component_image.SetGridValue(i, 1.0);
    }
  }
  
  // Find connected components (for salience and primitive markers)
  int *component_sizes = new int [ component_image.NEntries() ];
  int *component_memberships = new int [ component_image.NEntries() ];
  int ncomponents = component_image.ConnectedComponents(0.5, component_image.NEntries(),
    NULL, component_sizes, component_memberships);

  // Create features
  FETDescriptor descriptor;
  for (int i = 0; i < depth_image.XResolution(); i++) {
    for (int j = 0; j < depth_image.YResolution(); j++) {
      FETFeature *feature = NULL;
      int load_every_kth_pixel = 1;
      for (int s = 0; (s < load_every_kth_pixel) && !feature; s++) {
        for (int t = 0; (t < load_every_kth_pixel) && !feature; t++) {
          // Get index of pixel 
          int grid_index;
          int ix = i + s; if (ix >= depth_image.XResolution()) break;
          int iy = j + t; if (iy >= depth_image.YResolution()) break;
          depth_image.IndicesToIndex(ix, iy, grid_index);

          // Check boundary
          int boundary = (int) (boundary_image.GridValue(grid_index) + 0.5);
          if (boundary == 0) continue;

          // Check depth
          RNScalar depth = depth_image.GridValue(grid_index);
          if (RNIsNegativeOrZero(depth)) continue;
          if ((max_depth > 0) && (depth > max_depth)) continue;

          // Compute generator
          int generator_type = UNKNOWN_FEATURE_TYPE;
          if (boundary == R3_SURFEL_BORDER_BOUNDARY_FLAG) generator_type = BORDER_FEATURE_TYPE;
          else if (boundary == R3_SURFEL_SILHOUETTE_BOUNDARY_FLAG) generator_type = SILHOUETTE_FEATURE_TYPE;
          else if (boundary == R3_SURFEL_SHADOW_BOUNDARY_FLAG) generator_type = SHADOW_FEATURE_TYPE;
          else RNAbort("Unknown boundary type: $d\n", boundary);
          
          // Compute shape type
          int shape_type = PLANE_FEATURE_SHAPE;
          if (generator_type == SILHOUETTE_FEATURE_TYPE) shape_type = LINE_FEATURE_SHAPE;
          else if (generator_type == SHADOW_FEATURE_TYPE) shape_type = LINE_FEATURE_SHAPE;

          // Compute salience
          RNScalar salience = 1.0;
          if (generator_type == BORDER_FEATURE_TYPE) salience = 1E-3;
          else if (generator_type == SHADOW_FEATURE_TYPE) salience = 1E-3;
          else if (generator_type == SILHOUETTE_FEATURE_TYPE) {
            if (ncomponents > 0) {
              int component = component_memberships[grid_index];
              int component_size = component_sizes[component];
              RNScalar target_component_size = 32.0;
              salience = component_size / target_component_size;            
              if (salience > 1.0) salience = 1.0;
              if (depth > 1) salience /= depth;
            }
          }

          // Create feature
          feature = CreateFeature(reconstruction, shape, ix, iy, shape_type, generator_type, descriptor,
            px_image, py_image, pz_image, nx_image, ny_image, nz_image,
            depth_image, radius_image, boundary_image, color_image, salience,
            viewpoint, towards, up, previous_features, update_previous_features);

          // Set primitive marker
          if (feature) feature->primitive_marker = component_memberships[grid_index];
        }
      }
    }
  }

  // Delete info about silhouette image connected components
  delete [] component_sizes;
  delete [] component_memberships;

  // Return success
  return 1;
}



static int
CreateUniformFeatures(FETReconstruction *reconstruction, FETShape *shape,
  const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image, 
  const R2Grid& boundary_image, const R2Image& color_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  R3Kdtree<FETFeature *>& previous_features, RNBoolean update_previous_features = TRUE)
{
  // Compute salience
  RNScalar salience = 0.01;
  
  // Create features
  FETDescriptor descriptor;
  for (int i = 0; i < depth_image.XResolution(); i++) {
    for (int j = 0; j < depth_image.YResolution(); j++) {
      if (boundary_image.GridValue(i, j) != 0.0) continue;
      CreateFeature(reconstruction, shape, i, j, PLANE_FEATURE_SHAPE, UNIFORM_FEATURE_TYPE, descriptor,
        px_image, py_image, pz_image, nx_image, ny_image, nz_image,
        depth_image, radius_image, boundary_image, color_image, salience,
        viewpoint, towards, up, previous_features, update_previous_features);
    }
  }

  // Return success
  return 1;
}



static FETFeature *
CreatePlaneFeature(FETReconstruction *reconstruction, FETShape *shape, Segmentation *segmentation,
  Cluster *cluster, Point *point, int primitive_marker,
  R3Kdtree<FETFeature *>& previous_features, RNBoolean update_previous_features = TRUE)
{
  // Check spacing
  if (min_feature_spacing > 0) {
    FETFeature *feature = previous_features.FindAny(point->position, 0, min_feature_spacing);
    if (feature) { feature->primitive_marker = primitive_marker; return NULL; }
  }
  
  // Compute radius
  RNScalar radius = point->radius;
  if (min_feature_spacing > 0) {
    RNLength min_radius = 0.5 * min_feature_spacing;
    if (radius < min_radius) radius = min_radius;
  }

  // Compute salience
  RNScalar salience = (cluster) ? 0.1 : 0.01;
  if (point->depth > 1) salience /= point->depth;

  // Compute flags
  RNFlags flags = point->boundary;
  if (cluster) flags.Add(FEATURE_IS_PLANAR);

  // Compute generator type
  int generator_type = (cluster) ? PLANE_FEATURE_TYPE : UNIFORM_FEATURE_TYPE;

  // Create feature
  FETFeature *feature = new FETFeature(reconstruction, PLANE_PRIMITIVE_TYPE,
    point->position, R3zero_vector, point->normal, radius, FETDescriptor(), point->color, flags);

  // Set feature properties
  feature->SetGeneratorType(generator_type);
  feature->SetSalience(salience);
  feature->primitive_marker = primitive_marker;

  // Insert feature into shape
  shape->InsertFeature(feature);

  // Update previous features
  if ((min_feature_spacing > 0) && (update_previous_features)) {
    previous_features.InsertPoint(feature);
  }
  
  // Return feature
  return feature;
}



static int
CreatePlaneFeatures(FETReconstruction *reconstruction, FETShape *shape, Segmentation *segmentation,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  R3Kdtree<FETFeature *>& previous_features, RNBoolean update_previous_features = TRUE)
{
  // Create features for each cluster
  for (int i = 0; i < segmentation->clusters.NEntries(); i++) {
    Cluster *cluster = segmentation->clusters.Kth(i);

    // Create boundary features
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      Point *point = cluster->points.Kth(j);
      if (point->boundary == 0) continue;
      CreatePlaneFeature(reconstruction, shape, segmentation, cluster, point, i, previous_features, update_previous_features);
    }

    // Create non-boundary features
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      Point *point = cluster->points.Kth(j);
      if (point->boundary != 0) continue;
      CreatePlaneFeature(reconstruction, shape, segmentation, cluster, point, i, previous_features, update_previous_features);
    }
  }

  // Create plane features for points not in any cluster
  for (int i = 0; i < segmentation->points.NEntries(); i++) {
    Point *point = segmentation->points.Kth(i);
    if (point->cluster) continue;
    CreatePlaneFeature(reconstruction, shape, segmentation, NULL, point, -1, previous_features, update_previous_features);
  }    
  
  // Return success
  return 1;
}



static int
CreateCreaseFeatures(FETReconstruction *reconstruction, FETShape *shape, Segmentation *segmentation,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  R3Kdtree<FETFeature *>& previous_features, RNBoolean update_previous_features = TRUE)
{
  // Convenient variables
  RNArray<Pair *> pairs;
  R3Vector right = towards % up; right.Normalize();
  RNScalar max_plane_intersection_normal_dot = cos(min_plane_intersection_normal_angle);
  const int min_cluster_points_for_crease = 32;
  
  // Consider each cluster
  for (int i = 0; i < segmentation->clusters.NEntries(); i++) {
    Cluster *cluster0 = segmentation->clusters.Kth(i);
    if (cluster0->points.NEntries() < min_cluster_points_for_crease) continue;
    Primitive *primitive0 = &cluster0->primitive;
    if (primitive0->primitive_type != PLANE_PRIMITIVE_TYPE) continue;
    const R3Plane& plane0 = primitive0->plane;
    const R3Vector& normal0 = plane0.Normal();
    if (RNIsZero(normal0.Length())) continue;

    // Consider each point in the cluster
    for (int i0 = 0; i0 < cluster0->points.NEntries(); i0++) {
      Point *point0 = cluster0->points.Kth(i0);
      const R3Point& position0 = point0->position;

      // Check if another feature is within minimum spacing
      if (min_feature_spacing > 0) {
        if (previous_features.FindAny(position0, 0, min_feature_spacing)) { continue; }
      }

      // Consider each neighbor point
      for (int i1 = 0; i1 < point0->neighbors.NEntries(); i1++) {
        Point *point1 = point0->neighbors.Kth(i1);
        if (point1->grid_index <= point0->grid_index) continue;
        Cluster *cluster1 = point1->cluster;
        if (!cluster1) continue;
        if (cluster1 == cluster0) continue;
        if (cluster1->points.NEntries() < min_cluster_points_for_crease) continue;
        Primitive *primitive1 = &cluster1->primitive;
        if (primitive1->primitive_type != PLANE_PRIMITIVE_TYPE) continue;
        const R3Point& position1 = point1->position;
        const R3Plane& plane1 = primitive1->plane;
        R3Vector normal1 = plane1.Normal();
        RNScalar dot = normal1.Dot(normal0);
        if (fabs(dot) > max_plane_intersection_normal_dot) continue;
        if (RNIsZero(normal1.Length())) continue;
        
        // Compute intersection line
        R3Line line;
        if (!R3Intersects(R3Plane(position0, normal0), R3Plane(position1, normal1), &line)) continue;

        // Compute position
        R3Point position = 0.5 * (position0 + position1);
        position.Project(line);

        // Compute direction
        R3Vector direction = line.Vector();
        if (direction.Dot(right) < -0.1) direction.Flip();
        if (direction.Dot(up) < -0.1) direction.Flip();

        // Compute radius
        RNScalar radius = 0.5 * (point0->radius + point1->radius);
        if (min_feature_spacing > 0) {
          RNLength min_radius = 0.5 * min_feature_spacing;
          if (radius < min_radius) radius = min_radius;
        }

        // Compute salience
        RNScalar salience0 = 0.1 * cluster0->total_affinity / min_cluster_points;
        RNScalar salience1 = 0.1 * cluster1->total_affinity / min_cluster_points;
        RNScalar salience = (salience0 < salience1) ? salience0 : salience1;
        if (salience > 1.0) salience = 1.0;
        RNScalar depth0 = point0->depth;
        RNScalar depth1 = point1->depth;
        RNScalar depth = 0.5 * (depth0 + depth1);
        if (depth > 1) salience /= depth;
        
        // Compute generator type
        int generator_type = VALLEY_FEATURE_TYPE;
        R3Point c = 0.5 * (primitive0->centroid + primitive1->centroid);
        if (R3SignedDistance(plane0, c) < 0) generator_type = RIDGE_FEATURE_TYPE;
        else if (R3SignedDistance(plane1, c) < 0) generator_type = RIDGE_FEATURE_TYPE;

        // Compute flags
        RNFlags flags = point0->boundary | point1->boundary | FEATURE_IS_LINEAR;

        // Compute primitive marker
        int primitive_marker = FindPairIndex(pairs, cluster0, cluster1);
        if (primitive_marker < 0) {
          primitive_marker = pairs.NEntries();
          pairs.Insert(new Pair(cluster0, cluster1));
        }

        // Compute normal
        R3Vector normal = normal0 + normal1;
        normal.Normalize();
        
        // Create feature
        FETFeature *feature = new FETFeature(reconstruction, LINE_PRIMITIVE_TYPE,
          position, direction, normal, radius, FETDescriptor(), point0->color, flags);

        // Set feature properties
        feature->SetGeneratorType(generator_type);
        feature->SetSalience(salience);
        feature->primitive_marker = primitive_marker;

        // Insert feature into shape
        shape->InsertFeature(feature);

        // Update previous features
        if ((min_feature_spacing > 0) && (update_previous_features)) {
          previous_features.InsertPoint(feature);
        }
      }
    }
  }

  // Delete pairs
  for (int i = 0; i < pairs.NEntries(); i++) delete pairs[i];

  // Return success
  return 1;
}



static int
LoadFeatures(FETReconstruction *reconstruction,
  const char *scan_name, Segmentation *segmentation, const R2Grid& depth_image, 
  const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& radius_image, const R2Grid& boundary_image, const R2Image& color_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up)
{
  // Compute bounding box
  R3Box bbox = R3null_box;
  for (int i = 0; i < depth_image.NEntries(); i++) {
    RNScalar depth = depth_image.GridValue(i);
    if (depth <= 0) continue; 
    if ((max_depth > 0) && (depth > max_depth)) continue;
    RNScalar px = px_image.GridValue(i);
    RNScalar py = py_image.GridValue(i);
    RNScalar pz = pz_image.GridValue(i);
    R3Point position(px, py, pz);
    bbox.Union(position);
  }

  // Initialize kdtree of features (for min_feature_spacing)
  FETFeature tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
  R3Kdtree<FETFeature *> previous_features(bbox, position_offset);

  // Create shape
  FETShape *shape = CreateShape(reconstruction, scan_name, viewpoint, towards, up);
  if (!shape) return 0;

  // Create boundary features
  if (create_boundary_features) {
    if (!CreateBoundaryFeatures(reconstruction, shape, 
      px_image, py_image, pz_image, nx_image, ny_image, nz_image, 
      depth_image, radius_image, boundary_image, color_image, 
      viewpoint, towards, up, previous_features)) return 0;
  }

  // Create crease features
  if (segmentation && create_crease_features) {
    if (!CreateCreaseFeatures(reconstruction, shape, segmentation,
      viewpoint, towards, up, previous_features)) return 0;
  }

  // Create plane features
  if (segmentation && create_plane_features) {
    if (!CreatePlaneFeatures(reconstruction, shape, segmentation,
      viewpoint, towards, up, previous_features)) return 0;
  }

  // Create uniform features
  if (create_uniform_features) {
    if (!CreateUniformFeatures(reconstruction, shape, 
      px_image, py_image, pz_image, nx_image, ny_image, nz_image, 
      depth_image, radius_image, boundary_image, color_image, 
      viewpoint, towards, up, previous_features)) return 0;
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Generic processing stuff
////////////////////////////////////////////////////////////////////////

static int
ProcessRGBDScan(R3SurfelScene *scene, FETReconstruction *reconstruction,
  const char *scan_name, const char *depth_name, const char *color_name,
  R3Matrix depth_intrinsics_matrix, R3Matrix color_intrinsics_matrix, const R4Matrix& extrinsics_matrix,
  const char *dataset_format = NULL)
{
  // Determine camera coordinate system
  R3Point viewpoint = extrinsics_matrix * R3zero_point;
  R3Vector towards = extrinsics_matrix * R3negz_vector;
  R3Vector up = extrinsics_matrix * R3posy_vector;
  
   // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("  Processing %s\n", scan_name);
    fflush(stdout);
  }
  
  // Print timing message
  RNTime step_time;
  if (print_debug) {
    printf("    A\n");
    fflush(stdout);
    step_time.Read();
  }

  // Read depth image
  R2Grid depth_image;
  char input_filename[4096];
  if (input_depth_directory) sprintf(input_filename, "%s/%s", input_depth_directory, depth_name);
  else sprintf(input_filename, "%s", depth_name);
  if (!ReadDepthImage(depth_image, input_filename, dataset_format)) return 0;

  // Print timing message
  if (print_debug) {
    printf("    B %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }

  // Get color image
  R2Image color_image;
  if (color_name && strcmp(color_name, "-")) {
    if (input_color_directory) sprintf(input_filename, "%s/%s", input_color_directory, color_name);
    else sprintf(input_filename, "%s", color_name);
    if (!ReadColorImage(color_image, input_filename, dataset_format)) return 0;
  }

  // Print timing message
  if (print_debug) {
    printf("    C %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }

  // Resample images
  if ((max_image_resolution > 0) && (max_image_resolution != depth_image.XResolution())) {
    int xresolution = max_image_resolution;
    int yresolution = (int) ((depth_image.YResolution() * xresolution / (double) depth_image.XResolution()) + 0.5);
    if (!ResampleDepthImage(depth_image, depth_intrinsics_matrix, xresolution, yresolution)) return 0;
  }
  if ((color_image.Width() != depth_image.XResolution()) || (color_image.Height() != depth_image.YResolution())) {
    if (!ResampleColorImage(color_image, color_intrinsics_matrix, depth_image.XResolution(), depth_image.YResolution())) return 0;
  }
  
  // Print timing message
  if (print_debug) {
    printf("    D %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }

  // Create boundary image
  R2Grid boundary_image;
  if (!CreateBoundaryImage(depth_image, boundary_image)) return 0;

  // Print timing message
  if (print_debug) {
    printf("    E %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }

  // Create position images
  R2Grid px_image, py_image, pz_image;
  if (!CreatePositionImages(depth_image, depth_intrinsics_matrix, extrinsics_matrix, px_image, py_image, pz_image)) return 0;

  // Print timing message
  if (print_debug) {
    printf("    F %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }

  // Create normal and radius images
  R2Grid nx_image, ny_image, nz_image, radius_image;
  if (!CreateNormalImages(depth_image, 
    px_image, py_image, pz_image, boundary_image, 
    viewpoint, towards, up,
    nx_image, ny_image, nz_image, radius_image)) return 0;
  
  // Print timing message
  if (print_debug) {
    printf("    G %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }

  // Create segmentation
  Segmentation *segmentation = NULL;
  if (create_plane_features || create_crease_features) {
    segmentation = CreateSegmentation(
      px_image, py_image, pz_image, nx_image, ny_image, nz_image, 
      depth_image, radius_image, boundary_image, color_image, 
      viewpoint, towards, up);
  }

  // Print timing message
  if (print_debug) {
    printf("    H %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }

  // Load surfels
  if (scene) {
    if (segmentation) {
      // Load surfels with segmentation
      LoadSurfels(scene, scan_name, segmentation, depth_image, 
        viewpoint, towards, up);
    }
    else {
      // Load surfels without segmentation
      LoadSurfels(scene, scan_name, segmentation, depth_image, 
        px_image, py_image, pz_image, nx_image, ny_image, nz_image, 
        radius_image, boundary_image, color_image, 
        viewpoint, towards, up);
    }
  }

  // Print timing message
  if (print_debug) {
    printf("    I %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }
  
  // Load reconstruction features
  if (reconstruction) {
    LoadFeatures(reconstruction, scan_name, segmentation, depth_image, 
      px_image, py_image, pz_image, nx_image, ny_image, nz_image, 
      radius_image, boundary_image, color_image, 
      viewpoint, towards, up);
  }

  // Print timing message
  if (print_debug) {
    printf("    J %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }

  // Write mesh
  if (output_mesh_directory) {
    char output_filename[1024];
    sprintf(output_filename, "mkdir -p %s", output_mesh_directory);
    system(output_filename);
    sprintf(output_filename, "%s/%s.ply", output_mesh_directory, scan_name);
    WriteMesh(output_filename, segmentation, depth_image, 
              px_image, py_image, pz_image, nx_image, ny_image, nz_image, 
              radius_image, boundary_image, color_image, 
              viewpoint, towards, up);
  }

  // Print timing message
  if (print_debug) {
    printf("    I %g\n", step_time.Elapsed());
    fflush(stdout);
    step_time.Read();
  }
  
  // Delete segmentation
  if (segmentation) delete segmentation;
  
  // Print statistics
  if (print_debug) {
    printf("    Total = %.2f\n", start_time.Elapsed());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ProcessRGBDConfigurationFile(const char *input_configuration_name,
  R3SurfelScene *scene, FETReconstruction *reconstruction)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  int load_count = 0;
  if (print_verbose) {
    printf("Processing images from %s ...\n", input_configuration_name);
    fflush(stdout);
  }

  // Open input configuration file
  FILE *input_configuration_fp = fopen(input_configuration_name, "r");
  if (!input_configuration_fp) {
    fprintf(stderr, "Unable to open configuration file %s\n", input_configuration_name);
    return 0;
  }

  // Initialize dataset_format and intrinsics matrix
  const char *dataset_format = "";
  R3Matrix depth_intrinsics_matrix = R3identity_matrix;
  R3Matrix color_intrinsics_matrix = R3identity_matrix;

  // Parse file
  char buffer[4096];
  int line_number = 0;
  while (fgets(buffer, 4096, input_configuration_fp)) {
    char cmd[4096];
    line_number++;
    if (sscanf(buffer, "%s", cmd) != (unsigned int) 1) continue;
    if (cmd[0] == '#') continue;

    // Check cmd
    if (!strcmp(cmd, "dataset")) {
      // Parse dataset format
      char dataset_name[4096];
      if (sscanf(buffer, "%s%s", cmd, dataset_name) != (unsigned int) 2) {
        fprintf(stderr, "Error parsing line %d of %s\n", line_number, input_configuration_name);
        return 0;
      }

      // Remember dataset format
      // sunrgbd, matterport, nyu, etc.
      char *dataset_namep = dataset_name;
      while (*dataset_namep) { *dataset_namep = tolower(*dataset_namep); dataset_namep++; }
      dataset_format = strdup(dataset_name);
    }
    else if (!strcmp(cmd, "intrinsics") || !strcmp(cmd, "depth_intrinsics") || !strcmp(cmd, "color_intrinsics")) {
      // Parse filename
      char filename[4096];
      if (sscanf(buffer, "%s%s", cmd, filename) != (unsigned int) 2) {
        fprintf(stderr, "Error parsing line %d of %s\n", line_number, input_configuration_name);
        return 0;
      }

      // Read intrinsics matrix file
      R3Matrix intrinsics_matrix;
      if (!ReadIntrinsicsMatrix(intrinsics_matrix, filename, dataset_format)) {
        fprintf(stderr, "Unable to read intrinsics file %s\n", filename); 
        return 0;
      }

      // Assig depth and/or depth intrinsics matrices
      if (strcmp(cmd, "color_intrinsics")) depth_intrinsics_matrix = intrinsics_matrix;
      if (strcmp(cmd, "depth_intrinsics")) color_intrinsics_matrix = intrinsics_matrix;
    }
    else if (!strcmp(cmd, "intrinsics_matrix") || !strcmp(cmd, "depth_intrinsics_matrix") || !strcmp(cmd, "color_intrinsics_matrix")) {
      // Parse matrix
      double m[9];
      if (sscanf(buffer, "%s%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], &m[8]) != (unsigned int) 10) {
        fprintf(stderr, "Error parsing line %d of %s\n", line_number, input_configuration_name);
        return 0;
      }

      // Assign both depth and color intrinsics matrix 
      if (strcmp(cmd, "color_intrinsics_matrix")) depth_intrinsics_matrix = R3Matrix(m);
      if (strcmp(cmd, "depth_intrinsics_matrix")) color_intrinsics_matrix = R3Matrix(m);
    }
    else if (!strcmp(cmd, "depth_directory")) {
      // Parse directory name
      char dirname[4096];
      if (sscanf(buffer, "%s%s", cmd, dirname) != (unsigned int) 2) {
        fprintf(stderr, "Error parsing line %d of %s\n", line_number, input_configuration_name);
        return 0;
      }

      // Assign directory name
      input_depth_directory = strdup(dirname);
    }
    else if (!strcmp(cmd, "color_directory") || !strcmp(cmd, "image_directory")) {
      // Parse directory name
      char dirname[4096];
      if (sscanf(buffer, "%s%s", cmd, dirname) != (unsigned int) 2) {
        fprintf(stderr, "Error parsing line %d of %s\n", line_number, input_configuration_name);
        return 0;
      }

      // Assign directory name
      input_color_directory = strdup(dirname);
    }
    else if (!strcmp(cmd, "scan")) {
      // Check if should load image
      static int scan_count = -1; scan_count++;
      if ((load_every_kth_image > 1) && ((scan_count % load_every_kth_image) != 0)) continue; 
      if (scan_count < load_images_starting_at_index) continue;
      if (scan_count > load_images_ending_at_index) continue;

      // Parse image name and alignment transformation
      RNScalar m[16];
      char depth_name[4096], color_name[4096];
      if (sscanf(buffer, "%s%s%s%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", cmd, 
         depth_name, color_name,
         &m[0], &m[1], &m[2], &m[3], &m[4], &m[5], &m[6], &m[7], 
         &m[8], &m[9], &m[10], &m[11], &m[12], &m[13], &m[14], &m[15]) != (unsigned int) 19) {
        fprintf(stderr, "Error parsing line %d of %s\n", line_number, input_configuration_name);
        return 0;
      }

      // Get scan name
      char scan_name_buffer[4096];
      strncpy(scan_name_buffer, depth_name, 4096);
      char *scan_name = strrchr(scan_name_buffer, '/');
      if (!scan_name) scan_name = scan_name_buffer;
      char *endp = strrchr(scan_name, '.');
      if (endp) *endp = '\0';

      // Compute extrinsics matrix
      R4Matrix extrinsics_matrix(m);
      if (!strcmp(dataset_format, "nyu")) {
        R4Matrix rotation_matrix(-1, 0, 0, 0,   0, 0, 1, 0,   0, 1, 0, 0,   0, 0, 0, 1);
        extrinsics_matrix = rotation_matrix * extrinsics_matrix;
      }

      // Read image
      if (!ProcessRGBDScan(scene, reconstruction, scan_name, depth_name, color_name, 
        depth_intrinsics_matrix, color_intrinsics_matrix, extrinsics_matrix,
        dataset_format)) return 0;

      // Update scan counter
      load_count++;
    }
  }

  // Close configuration file
  if (input_configuration_fp) fclose(input_configuration_fp);

  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Scans = %d\n", load_count);
    if (scene) {
      printf("  # Surfel Nodes = %d\n", scene->Tree()->NNodes());
      printf("  # Surfel Blocks = %d\n", scene->Tree()->Database()->NBlocks());
      printf("  # Surfels = %d\n", scene->Tree()->Database()->NSurfels());
    }
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Local variable to track parameter setting
  RNBoolean create_features = FALSE;

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-debug")) print_debug = 1;
      else if (!strcmp(*argv, "-progress")) print_progress = 1;
      else if (!strcmp(*argv, "-output_mesh_directory")) { argc--; argv++; output_mesh_directory =*argv; }
      else if (!strcmp(*argv, "-load_image_at_index")) { argc--; argv++; load_images_starting_at_index = load_images_ending_at_index = atoi(*argv); }
      else if (!strcmp(*argv, "-load_images_starting_at_index")) { argc--; argv++; load_images_starting_at_index = atoi(*argv); }
      else if (!strcmp(*argv, "-load_images_ending_at_index")) { argc--; argv++; load_images_ending_at_index = atoi(*argv); }
      else if (!strcmp(*argv, "-load_every_kth_image")) { argc--; argv++; load_every_kth_image = atoi(*argv); }
      else if (!strcmp(*argv, "-create_boundary_features")) { create_features = create_boundary_features = TRUE; }
      else if (!strcmp(*argv, "-create_uniform_features")) { create_features = create_uniform_features = TRUE; }
      else if (!strcmp(*argv, "-create_plane_features")) { create_features = create_plane_features = TRUE; }
      else if (!strcmp(*argv, "-create_crease_features")) { create_features = create_crease_features = TRUE; }
      else if (!strcmp(*argv, "-create_multiresolution_surfels")) create_multiresolution_surfels = TRUE;
      else if (!strcmp(*argv, "-normal_neighborhood_pixel_radius")) { argc--; argv++; normal_neighborhood_pixel_radius = atoi(*argv); }
      else if (!strcmp(*argv, "-normal_neighborhood_world_radius")) { argc--; argv++; normal_neighborhood_world_radius = atof(*argv); }
      else if (!strcmp(*argv, "-max_image_resolution")) { argc--; argv++; max_image_resolution = atoi(*argv); }
      else if (!strcmp(*argv, "-min_feature_spacing")) { argc--; argv++; min_feature_spacing = atof(*argv); }
      else if (!strcmp(*argv, "-max_depth")) { argc--; argv++; max_depth = atof(*argv); }
      else if (!strcmp(*argv, "-omit_corners")) omit_corners = 1;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_configuration_name) input_configuration_name = *argv;
      else if (!output_ssa_name && strstr(*argv, ".ssa")) output_ssa_name = *argv;
      else if (!output_ssb_name && strstr(*argv, ".ssb")) output_ssb_name = *argv;
      else if (!output_reconstruction_name && strstr(*argv, ".fca")) output_reconstruction_name = *argv;
      else if (!output_reconstruction_name && strstr(*argv, ".fcb")) output_reconstruction_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check filenames
  if (!input_configuration_name) {
    fprintf(stderr, "Usage: rgbdprocess inputconfigurationfile [outputdirectory] [ssafile ssbfile] [fetfile] [options]\n");
    return 0;
  }

  // Do everything if no program arguments
  if (!output_ssa_name && !output_ssb_name && !output_reconstruction_name) {
    // Get scene name
    char scene_name[4096];
    const char *strp = strrchr(input_configuration_name, '/');
    if (!strp) strp = input_configuration_name;
    strcpy(scene_name, strp);
    char *endp = strrchr(scene_name, '.');
    if (endp) *endp = '\0';

    // Get default output filenames
    char buffer[4096];
    sprintf(buffer, "%s.ssa", scene_name);
    output_ssa_name = strdup(buffer);
    sprintf(buffer, "%s.ssb", scene_name);
    output_ssb_name = strdup(buffer);
    sprintf(buffer, "%s.fcb", scene_name);
    output_reconstruction_name = strdup(buffer);
  }

  // Fill in default feature creation parameters
  if (output_reconstruction_name && !create_features) {
    create_boundary_features = TRUE;
    create_uniform_features = TRUE;
    create_plane_features = TRUE;
    create_crease_features = TRUE;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Open surfel scene
  R3SurfelScene *scene = NULL;
  if (output_ssa_name && output_ssb_name) {
    scene = OpenSurfelScene(output_ssa_name, output_ssb_name);
    if (!scene) exit(-1);
  }

  // Create reconstruction
  FETReconstruction *reconstruction = NULL;
  if (output_reconstruction_name) {
    reconstruction = CreateReconstruction();
    if (!reconstruction) exit(-1);
  }

  // Process RGBD data
  if (!ProcessRGBDConfigurationFile(input_configuration_name, scene, reconstruction)) exit(-1);

  // Write reconstruction
  if (reconstruction && output_reconstruction_name) {
    if (!WriteReconstruction(reconstruction, output_reconstruction_name)) exit(-1);
  }

  // Process surfel scene
  if (scene && create_multiresolution_surfels) {
    if (!CreateMultiresolutionSurfels(scene)) exit(-1);
  }

  // Close surfel scene
  if (scene && output_ssa_name && output_ssb_name) {
    if (!CloseSurfelScene(scene)) exit(-1);
  }

  // Return success 
  return 0;
}



