////////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////

static int min_cluster_points = 10;
static int min_clusters = 0;
static int max_clusters = 0;
static double min_cluster_coverage = 0;
static double max_cluster_diameter = 16;
static double max_cluster_primitive_distance = 0.1;
static double max_cluster_normal_angle = RN_PI / 4.0;
static double max_neighbor_distance_factor = 16;
static double max_pair_centroid_distance = 16;
static double max_pair_primitive_distance = 0.1;
static double max_pair_normal_angle = RN_PI / 4.0;
static double min_pair_affinity = 1.0E-6;
static double min_plane_intersection_normal_angle = RN_PI/4.0;
// static int max_refinement_iterations = 2;  // could be 1 with little change in quality
// static int max_ransac_iterations = 1; // was 4, which may be better for SUN3D
static int initialize_hierarchically = TRUE; 
static int print_progress = FALSE;

static int max_refinement_iterations = 1;  // WORKS FINE
static int max_ransac_iterations = 0;



///////////////////////////////////////////////////////////////////////
// Shape types
////////////////////////////////////////////////////////////////////////

enum {
  NULL_PRIMITIVE_TYPE,
  POINT_PRIMITIVE_TYPE,
  LINE_PRIMITIVE_TYPE,
  PLANE_PRIMITIVE_TYPE,
  PLANAR_GRID_PRIMITIVE_TYPE,
  NUM_PRIMITIVE_TYPES
};



////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////

struct Point {
public:
  Point(void);
public:
  RNScalar depth;
  R3Point position;
  R3Vector normal;
  RNScalar radius;
  RNRgb color;
  unsigned int boundary;
  RNArray<Point *> neighbors;
  struct Cluster *cluster;
  RNScalar cluster_affinity;
  int cluster_index;
  int grid_index;
  int mark;
};

struct Primitive {
  Primitive(int primitive_type = 0);
  Primitive(const Primitive& primitive);
  Primitive(Point *seed_point, const RNArray<Point *> *points = NULL);
  RNLength Distance(const R3Point& position) const;
  void Update(const R3Point& point);
  void Update(const R3Line& line);
  void Update(const R3Plane& plane);
  void Update(Point *seed_point = NULL, const RNArray<Point *> *points = NULL);
  void Update(Primitive primitive1, Primitive primitive2, RNScalar weight1 = 1.0, RNScalar weight2 = 1.0);
public:
  int primitive_type;
  R3Box bbox;
  R3Point centroid;
  R3Line line;
  R3Plane plane;
};

struct Cluster {
public:
  Cluster(Point *seed_point = NULL, int primitive_type = 0);
  Cluster(Point *seed_point, const Primitive& primitive);
  Cluster(Cluster *child1, Cluster *child2);
  ~Cluster(void);
  RNScalar Coverage(void);
  void EmptyPoints(void);
  void InsertPoint(Point *point, RNScalar affinity = 1.0);
  void RemovePoint(Point *point);
  void InsertChild(Cluster *child);
  void RemoveChild(Cluster *child);
  int UpdatePoints(const R3Kdtree<Point *> *kdtree);
  int UpdatePrimitive(void);
  RNScalar Affinity(Point *point) const;
  RNScalar Affinity(Cluster *cluster) const;
public:
  Point *seed_point;
  RNArray<Point *> points;
  Cluster *parent;
  RNArray<Cluster *> children;
  RNArray<struct Pair *> pairs;
  Primitive primitive;
  RNScalar possible_affinity; 
  RNScalar total_affinity;
  int segmentation_index;
};

struct Pair {
public:
  Pair(Cluster *cluster1 = NULL, Cluster *cluster2 = NULL, RNScalar affinity = 0);
  ~Pair(void);
public:
  Cluster *clusters[2];
  int cluster_index[2];
  RNScalar affinity; 
  Pair **heapentry;
};

struct Segmentation {
public:
  Segmentation(void);
  ~Segmentation(void);
  RNScalar Affinity(void) const;
  int NUnclusteredPoints(void) const;
public:
  int CreatePoints(const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image, 
    const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
    const R2Grid& depth_image, const R2Grid& radius_image,
    const R2Grid& boundary_image, const R2Image& color_image);
  int CreateClusters(int primitive_type);
  int CreateSingletonClusters(int primitive_type);
  int CreateRansacClusters(int primitive_type);
  int RefineClusters(void);  
  int DeleteClusters(void);  
  int MergeClusters(void);  
  int SplitClusters(void);
public:
  int ReadClusterImage(const char *filename);
  int WriteClusterImage(int xres, int yres, const char *filename) const;
public:
  RNArray<Point *> points;
  R3Kdtree<Point *> *kdtree;
  RNArray<Cluster *> clusters;
  Point *point_buffer;
};




////////////////////////////////////////////////////////////////////////
// Point member functions
////////////////////////////////////////////////////////////////////////

Point::
Point(void)
  : depth(-1),
    position(0,0,0),
    normal(0,0,0),
    radius(0),
    boundary(0),
    neighbors(),
    cluster(NULL),
    cluster_affinity(0),
    cluster_index(-1),
    grid_index(-1),
    mark(0)
{
}



////////////////////////////////////////////////////////////////////////
// Primitive member functions
////////////////////////////////////////////////////////////////////////

Primitive::
Primitive(int primitive_type)
  : primitive_type(primitive_type),
    bbox(R3null_box),
    centroid(R3zero_point),
    line(R3null_line),
    plane(R3null_plane)
{
}



Primitive::
Primitive(const Primitive& primitive)
  : primitive_type(primitive.primitive_type),
    bbox(primitive.bbox),
    centroid(primitive.centroid),
    line(primitive.line),
    plane(primitive.plane)
{
}



Primitive::
Primitive(Point *seed_point, const RNArray<Point *> *points)
  : primitive_type(NULL_PRIMITIVE_TYPE),
    bbox(R3null_box),
    centroid(R3zero_point),
    line(R3null_line),
    plane(R3null_plane)
{
  // Initialize primitive based on points
  Update(seed_point, points); 
}



RNLength Primitive::
Distance(const R3Point& position) const
{
  // Return distance from primitive to point
  if (primitive_type == POINT_PRIMITIVE_TYPE) return R3Distance(centroid, position);
  else if (primitive_type == LINE_PRIMITIVE_TYPE) return R3Distance(line, position);
  else if (primitive_type == PLANE_PRIMITIVE_TYPE) return R3Distance(plane, position);
#if 0
  else if (primitive_type == PLANAR_GRID_PRIMITIVE_TYPE) {
    R2Point grid_position = planar_grid.GridPosition(position);
    int ix = (int) (grid_position.X() + 0.5);
    if ((ix < 0) || (ix >= planar_grid.XResolution())) return RN_INFINITY;
    int iy = (int) (grid_position.Y() + 0.5);
    if ((iy < 0) || (iy >= planar_grid.YResolution())) return RN_INFINITY;
    RNScalar value = planar_grid.GridValue(ix, iy);
    if (value <= 0) return RN_INFINITY;
    return R3Distance(planar_grid.Plane(), position);
  }
#endif
  else {
    RNAbort("Unrecognized primitive type");
    return RN_INFINITY;
  }
}



void Primitive::
Update(const R3Point& point)
{
  // Set everything
  primitive_type = POINT_PRIMITIVE_TYPE;
  this->centroid = point;
  line = R3null_line;
  plane = R3null_plane;
}



void Primitive::
Update(const R3Line& line)
{
  // Set everything
  primitive_type = LINE_PRIMITIVE_TYPE;
  centroid = R3zero_point;
  centroid.Project(line);
  this->line = line;
  plane = R3null_plane;
}



void Primitive::
Update(const R3Plane& plane)
{
  // Set everything
  primitive_type = PLANE_PRIMITIVE_TYPE;
  centroid = R3zero_point;
  centroid.Project(plane);
  line = R3null_line;
  this->plane = plane;
}



void Primitive::
Update(Point *seed_point, const RNArray<Point *> *points)
{
  // Remember stuff about primitive (so can set same orientation)
  R3Vector previous_vector(0,0,0);
  if (primitive_type == LINE_PRIMITIVE_TYPE) previous_vector = line.Vector();
  if (primitive_type == PLANE_PRIMITIVE_TYPE) previous_vector = plane.Normal();

  // Update bounding box
  bbox = R3null_box;
  if (points) {
    for (int i = 0; i < points->NEntries(); i++) {
      Point *point = points->Kth(i);
      bbox.Union(point->position);
    }
  }

  // Initialize everything
  if (seed_point) {
    R3Point seed_position = seed_point->position;
    centroid = seed_position;
    line.Reset(seed_position, line.Vector()); 
    plane.Reset(seed_position, seed_point->normal);
    bbox.Union(seed_position);
  }
  else {
    // Temporary
    RNAbort("Need seed point");
  }

  // Update based on points
  if (points && (points->NEntries() > 0)) {
    // Allocate arrays of point positions and weights
    const int max_positions = 1024;
    R3Point *positions = new R3Point [ max_positions ];
    RNScalar *weights = new RNScalar [ max_positions ];

    // Fill arrays of point positions and weights
    int npositions = 0;
    int skip = points->NEntries() / max_positions + 1;
    for (int i = 0; i < points->NEntries(); i += skip) {
      Point *point = points->Kth(i);
      if (npositions >= max_positions-1) break;
      positions[npositions] = point->position;
      if (primitive_type == PLANE_PRIMITIVE_TYPE) weights[npositions] = fabs(plane.Normal().Dot(point->normal));
      else if (primitive_type == LINE_PRIMITIVE_TYPE) weights[npositions] = 1.0 - fabs(plane.Normal().Dot(point->normal));
      else weights[npositions] = 1.0;
      npositions++;
    }

    // Add seed point with 20% of the total weight
    if (seed_point) {
      positions[npositions] = seed_point->position;
      weights[npositions] = 0.2 * points->NEntries();
      npositions++;
    }

    // Compute centroid
    centroid = R3Centroid(npositions, positions, weights);

    // Update primitive parameters
    if ((primitive_type == NULL_PRIMITIVE_TYPE) && (npositions >= 2)) {
      RNScalar variances[3];
      R3Triad axes = R3PrincipleAxes(centroid, npositions, positions, weights, variances);
      if (variances[0] > RN_EPSILON) {
        if (variances[1] > RN_EPSILON) {
          RNScalar ratio10 = variances[1] / variances[0];
          RNScalar ratio21 = variances[2] / variances[1];
          if (ratio10 < ratio21) {
            primitive_type = LINE_PRIMITIVE_TYPE;
            line.Reset(centroid, axes[0]);
          }
          else {
            primitive_type = PLANE_PRIMITIVE_TYPE;
            plane.Reset(centroid, axes[2]);
          }
        }
      }
    }
    else if ((primitive_type == LINE_PRIMITIVE_TYPE) && (npositions >= 2)) {
      // Compute principle directions
      RNScalar variances[3];
      R3Triad axes = R3PrincipleAxes(centroid, npositions, positions, weights, variances);
      if (variances[0] > RN_EPSILON) {
        // Update line
        R3Vector direction = axes[0];
        line.Reset(centroid, direction);

        // Check if should flip line
        RNScalar dot = direction.Dot(previous_vector);
        if (dot < 0) line.Flip();
      }
    }
    else if ((primitive_type == PLANE_PRIMITIVE_TYPE) && (npositions >= 3)) {
      // Compute principle directions
      RNScalar variances[3];
      R3Triad axes = R3PrincipleAxes(centroid, npositions, positions, weights, variances);
      if (variances[1] > RN_EPSILON) {
        // Update plane
        R3Vector normal = axes[2];
        plane.Reset(centroid, normal);

        // Check if should flip plane
        if (seed_point) {
          RNScalar dot = normal.Dot(seed_point->normal);
          if (dot < 0) plane.Flip();
        }
        else {
          RNScalar dot = normal.Dot(previous_vector);
          if (dot < 0) plane.Flip();
        }
      }
    }

#if 0
    // Rasterize planar grid
    if (primitive_type == PLANAR_GRID_PRIMITIVE_TYPE) {
      // Rasterize points into planar grid
      R3PlanarGrid density(plane, bbox, min_cluster_spacing);
      for (int i = 0; i < points->NEntries(); i++) {
        Point *point = points->Kth(i);
        R3Point position = point->position;
        density.RasterizeWorldPoint(position.X(), position.Y(), position.Z(), 1.0);
      }

      // Reset planar grid
      planar_grid.Reset(plane, bbox, min_cluster_spacing);
      R3Point seed_position = seed_point->position;
      R2Point grid_position = planar_grid.GridPosition(seed_position);
      int ix = (int) (grid_position.X() + 0.5);
      int iy = (int) (grid_position.Y() + 0.5);
      FloodCopy(density.grid, planar_grid.grid, ix, iy);
    }
#endif

    // Delete array of point positions and weights
    delete [] positions;
    delete [] weights;
  }
}



void Primitive::
Update(Primitive primitive1, Primitive primitive2, RNScalar weight1, RNScalar weight2)
{
  // Just checking
  if (weight1 == 0) {
    primitive_type = primitive2.primitive_type;
    bbox = primitive2.bbox;
    centroid = primitive2.centroid;
    line = primitive2.line;
    plane = primitive2.plane;
  }
  else if (weight2 == 0) {
    primitive_type = primitive1.primitive_type;
    bbox = primitive1.bbox;
    centroid = primitive1.centroid;
    line = primitive1.line;
    plane = primitive1.plane;
  }
  else {
    // Update primitive type
    if (primitive1.primitive_type > primitive2.primitive_type) {
      primitive_type = primitive1.primitive_type;
      weight2 = 0;
    }
    else if (primitive2.primitive_type > primitive1.primitive_type) {
      primitive_type = primitive2.primitive_type;
      weight1 = 0;
    }
    else {
      primitive_type = primitive1.primitive_type;
    }

    // Update centroid
    centroid = R3zero_point;
    centroid += weight1 * primitive1.centroid;
    centroid += weight2 * primitive2.centroid;
    centroid /= weight1 + weight2;

    // Update bbox
    bbox = R3null_box;
    bbox.Union(primitive1.bbox);
    bbox.Union(primitive2.bbox);

    // Update other stuff
    line = R3null_line;
    plane = R3null_plane;
    if (primitive_type == LINE_PRIMITIVE_TYPE) {
      // Compute line
      R3Vector vector1 = primitive1.line.Vector();
      R3Vector vector2 = primitive2.line.Vector();
      if (vector1.Dot(vector2) < 0) vector2.Flip();
      R3Vector vector = R3zero_vector;
      vector += weight1 * vector1;
      vector += weight2 * vector2;
      vector /= weight1 + weight2;
      vector.Normalize();
      line.Reset(centroid, vector);
    }
    else if (primitive_type == PLANE_PRIMITIVE_TYPE) {
      // Compute plane
      R3Vector normal1 = primitive1.plane.Normal();
      R3Vector normal2 = primitive2.plane.Normal();
      if (normal1.Dot(normal2) < 0) normal2.Flip();
      R3Vector normal = R3zero_vector;
      normal += weight1 * normal1;
      normal += weight2 * normal2;
      normal /= weight1 + weight2;
      normal.Normalize();
      plane.Reset(centroid, normal);
    }
  }
}



////////////////////////////////////////////////////////////////////////
// Cluster member functions
////////////////////////////////////////////////////////////////////////

Cluster::
Cluster(Point *seed_point, int primitive_type)
  : seed_point(seed_point),
    points(),
    parent(NULL),
    children(),
    pairs(),
    primitive(primitive_type),
    possible_affinity(0),
    total_affinity(0),
    segmentation_index(-1)
{
  // Update primitive
  if (seed_point) primitive.Update(seed_point);
}



Cluster::
Cluster(Point *seed_point, const Primitive& primitive)
  : seed_point(seed_point),
    points(),
    parent(NULL),
    children(),
    pairs(),
    primitive(primitive),
    possible_affinity(0),
    total_affinity(0),
    segmentation_index(-1)
{
}



Cluster::
Cluster(Cluster *child1, Cluster *child2)
  : seed_point(NULL),
    points(),
    parent(NULL),
    children(),
    pairs(),
    primitive(),
    possible_affinity(0),
    total_affinity(0),
    segmentation_index(-1)
{
  // Assign seed point
  seed_point = child1->seed_point;

  // Update primitive
  primitive.Update(child1->primitive, child2->primitive, child1->points.NEntries(), child2->points.NEntries());

  // Insert points from child1
  while (!child1->points.IsEmpty()) {
    Point *point = child1->points.Tail();
    child1->RemovePoint(point);
    RNScalar affinity = Affinity(point);
    // RNScalar affinity = point->cluster_affinity; // THIS IS WRONG, USING AFFINITY TO OLD CLUSTER FOR SPEED
    if (affinity < 0) affinity = 0;
    possible_affinity += affinity;
    InsertPoint(point, affinity);
  }

  // Insert points from child2
  while (!child2->points.IsEmpty()) {
    Point *point = child2->points.Tail();
    child2->RemovePoint(point);
    RNScalar affinity = Affinity(point);
    // RNScalar affinity = point->cluster_affinity; // THIS IS WRONG, USING AFFINITY TO OLD CLUSTER FOR SPEED
    if (affinity < 0) affinity = 0;
    possible_affinity += affinity;
    InsertPoint(point, affinity);
  }

  // Update hierarchy
  child1->parent = this;
  child2->parent = this;
  children.Insert(child1);
  children.Insert(child2);
}



Cluster::
~Cluster(void)
{
  // Delete children
  // for (int i = 0; i < children.NEntries(); i++) {
  //   delete children.Kth(i);
  // }

  // Remove from parent
  if (parent) parent->children.Remove(this);

  // Empty points
  EmptyPoints();
}



RNScalar Cluster::
Coverage(void)
{
  // Return metric of how well cluster covers points
  if (possible_affinity == 0) return 0;
  return total_affinity / possible_affinity;
}



void Cluster::
EmptyPoints(void)
{
  // Update points
  for (int i = 0; i < points.NEntries(); i++) {
    Point *point = points.Kth(i);
    point->cluster = NULL;
    point->cluster_affinity = 0;
    point->cluster_index = -1;
  }

  // Empty points
  points.Empty();

  // Update affinity
  total_affinity = 0;
}



void Cluster::
InsertPoint(Point *point, RNScalar affinity)
{
  // Remove from previous cluster
  if (point->cluster ) {
    if (point->cluster == this) return;
    else point->cluster->RemovePoint(point);
  }

  // Update point
  point->cluster = this;
  point->cluster_index = points.NEntries();
  point->cluster_affinity = affinity;

  // Insert point
  points.Insert(point);

  // Update cluster
  total_affinity += point->cluster_affinity;
}



void Cluster::
RemovePoint(Point *point)
{
  // Just checking
  assert(point->cluster == this);
  assert(point->cluster_index >= 0);

  // Update cluster
  total_affinity -= point->cluster_affinity;

  // Remove point
  RNArrayEntry *entry = points.KthEntry(point->cluster_index);
  Point *tail = points.Tail();
  tail->cluster_index = point->cluster_index;
  points.EntryContents(entry) = tail;
  points.RemoveTail();

  // Update point
  point->cluster = NULL;
  point->cluster_index = -1;
  point->cluster_affinity = 0;
}



void Cluster::
InsertChild(Cluster *child)
{
  // Update primitive
  primitive.Update(this->primitive, child->primitive, this->points.NEntries(), child->points.NEntries());

  // Update affinities for current points
  if (points.NEntries() < 4 * child->points.NEntries()) {
    for (int i = 0; i < points.NEntries(); i++) {
      Point *point = points.Kth(i);
      RNScalar affinity = Affinity(point);
      if (affinity < 0) affinity = 0;
      possible_affinity += affinity - point->cluster_affinity;
      point->cluster_affinity = affinity;
    }
  }

  // Insert points from child
  while (!child->points.IsEmpty()) {
    Point *point = child->points.Tail();
    child->RemovePoint(point);
    RNScalar affinity = Affinity(point);
    if (affinity < 0) affinity = 0;
    possible_affinity += affinity;
    InsertPoint(point, affinity);
  }

  // Update hierarchy
  child->parent = this;
  children.Insert(child);
}



void Cluster::
RemoveChild(Cluster *child)
{
  // Remove child
  this->children.Remove(child);
  child->parent = NULL;
}



int Cluster::
UpdatePoints(const R3Kdtree<Point *> *kdtree)
{
  // Empty points
  // If do this, some points may end up as part of no cluster
  // If don't do this, some clusters may end up with outlier points if primitive changes a lot
  // EmptyPoints();  

  // Find points near primitive
  RNArray<Point *> points1;
  if (seed_point) {
    // Find connected set of points near primitive
    static int mark = 1;
    RNArray<Point *> stack;
    stack.Insert(seed_point);
    seed_point->mark = ++mark;
    while (!stack.IsEmpty()) {
      Point *point = stack.Tail();
      stack.RemoveTail();
      points1.Insert(point);
      for (int i = 0; i < point->neighbors.NEntries(); i++) {
        Point *neighbor = point->neighbors.Kth(i);
        if (neighbor->mark == mark) continue;
        neighbor->mark = mark;
        RNLength d = primitive.Distance(neighbor->position);
        if (d > max_cluster_primitive_distance) continue;
        stack.Insert(neighbor);
      }
    }
  }
  else if (kdtree) {
    // Find all points near primitive
    if (primitive.primitive_type == POINT_PRIMITIVE_TYPE) kdtree->FindAll(primitive.centroid, 0, max_cluster_primitive_distance, points1);
    else if (primitive.primitive_type == LINE_PRIMITIVE_TYPE) kdtree->FindAll(primitive.line, 0, max_cluster_primitive_distance, points1);
    else if (primitive.primitive_type == PLANE_PRIMITIVE_TYPE) kdtree->FindAll(primitive.plane, 0, max_cluster_primitive_distance, points1);
    else RNAbort("Unrecognized primitive type");
  }

  // Check points
  if (points1.NEntries() < min_cluster_points) {
    return 0;
  }
  
  // Allocate affinities
  RNScalar *affinities = new RNScalar [ points1.NEntries() ];
  if (!affinities) return 0;

  // Compute affinities
  possible_affinity = 0;
  RNArray<Point *> points2;
  for (int i = 0; i < points1.NEntries(); i++) {
    Point *point = points1.Kth(i);
    RNScalar affinity = Affinity(point);
    if (affinity <= 0) continue;
    affinities[points2.NEntries()] = affinity;
    points2.Insert(point);
    possible_affinity += affinity;
  }

  // Check points
  if (points2.NEntries() < min_cluster_points) {
    delete [] affinities;
    return 0;
  }

  // Compute assignments
  int npoints = 0;
  for (int i = 0; i < points2.NEntries(); i++) {
    Point *point = points2.Kth(i);
    if ((point->cluster) && (point->cluster != this)) {
      if (point->cluster->possible_affinity > 0) {
        RNScalar factor = possible_affinity / point->cluster->possible_affinity;
        if (point->cluster_affinity > factor * factor * affinities[i]) continue;
      }
    }
    npoints++;
  }

  // Check assignments
  if (npoints < min_cluster_points) {
    delete [] affinities;
    return 0;
  }

  // Insert points (should match previous loop)
  for (int i = 0; i < points2.NEntries(); i++) {
    Point *point = points2.Kth(i);
    if ((point->cluster) && (point->cluster != this)) {
      if (point->cluster->possible_affinity > 0) {
        RNScalar factor = possible_affinity / point->cluster->possible_affinity;
        if (point->cluster_affinity > factor * factor * affinities[i]) continue;
      }
    }
    InsertPoint(point, affinities[i]);
  }

  // Delete affinities
  delete [] affinities;

  // Return success
  return 1;
}



int Cluster::
UpdatePrimitive(void)
{
  // Update primitive
  primitive.Update(seed_point, &points);
  if (primitive.primitive_type == NULL_PRIMITIVE_TYPE) return 0;
  else return 1;
}



RNScalar Cluster::
Affinity(Point *point) const
{
  // Initialize affinity
  RNScalar affinity = 1.0;

  // Get useful variables
  R3Point position = point->position;

  // Check primitive distance 
  if (max_cluster_primitive_distance > 0) {
    RNLength primitive_distance = primitive.Distance(position);
    if (primitive_distance > max_cluster_primitive_distance) return 0;
    RNScalar primitive_distance_affinity = 1.0 - primitive_distance / max_cluster_primitive_distance;
    affinity *= primitive_distance_affinity;
  }

  // Check centroid distance
  if (max_cluster_diameter > 0) {
    RNLength centroid_distance = R3Distance(primitive.centroid, position);
    if (centroid_distance > max_cluster_diameter) return 0;
    RNScalar centroid_distance_affinity = 1.0 - centroid_distance / max_cluster_diameter;
    affinity *= centroid_distance_affinity;
  }

  // Check normal angle
  if (max_cluster_normal_angle > 0) {
    if (primitive.primitive_type == LINE_PRIMITIVE_TYPE) {
      RNScalar dot = fabs(primitive.line.Vector().Dot(point->normal));
      RNAngle normal_angle = (dot < 1) ? RN_PI_OVER_TWO - acos(dot) : RN_PI_OVER_TWO;
      if (normal_angle > max_cluster_normal_angle) return 0;
      RNScalar normal_angle_affinity = 1.0 - normal_angle / max_cluster_normal_angle;
      affinity *= normal_angle_affinity;
    }
    else if (primitive.primitive_type == PLANE_PRIMITIVE_TYPE) {
      RNScalar dot = fabs(primitive.plane.Normal().Dot(point->normal));
      RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
      if (normal_angle > max_cluster_normal_angle) return 0;
      RNScalar normal_angle_affinity = 1.0 - normal_angle / max_cluster_normal_angle;
      affinity *= normal_angle_affinity;
    }
  }

  // Return affinity
  return affinity;
}



RNScalar Cluster::
Affinity(Cluster *cluster) const
{
  // Initialize affinity
  RNScalar affinity = 1;

  // Compute centroid distance
  if (max_pair_centroid_distance > 0) {
    RNLength cluster_distance = R3Distance(primitive.centroid, cluster->primitive.centroid);
    if (cluster_distance > max_pair_centroid_distance) return 0;
    RNScalar cluster_distance_factor = 1.0 - cluster_distance / max_pair_centroid_distance;
    affinity *= cluster_distance_factor;
  }

  // Compute primitive distances
  if (max_pair_primitive_distance > 0) {
    // Compute point0-primitive1 distance
    RNLength primitive0_distance = primitive.Distance(cluster->primitive.centroid);
    if (primitive0_distance > max_pair_primitive_distance) return 0;
    RNScalar primitive0_distance_factor = 1.0 - primitive0_distance / max_pair_primitive_distance;
    affinity *= primitive0_distance_factor;

    // Compute point1-primitive0 distance
    RNLength primitive1_distance = cluster->primitive.Distance(primitive.centroid);
    if (primitive1_distance > max_pair_primitive_distance) return 0;
    RNScalar primitive1_distance_factor = 1.0 - primitive1_distance / max_pair_primitive_distance;
    affinity *= primitive1_distance_factor;
  }

  // Compute normal angle
  if (max_pair_normal_angle > 0) {
    if ((primitive.primitive_type == LINE_PRIMITIVE_TYPE) && (cluster->primitive.primitive_type == LINE_PRIMITIVE_TYPE)) {
      RNScalar dot = fabs(primitive.line.Vector().Dot(cluster->primitive.line.Vector()));
      RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
      if (normal_angle > max_cluster_normal_angle) return 0;
      RNScalar normal_angle_affinity = 1.0 - normal_angle / max_cluster_normal_angle;
      affinity *= normal_angle_affinity;
    }
    else if ((primitive.primitive_type == PLANE_PRIMITIVE_TYPE) && (cluster->primitive.primitive_type == LINE_PRIMITIVE_TYPE)) {
      RNScalar dot = fabs(primitive.plane.Normal().Dot(cluster->primitive.line.Vector()));
      RNAngle normal_angle = (dot < 1) ? RN_PI_OVER_TWO - acos(dot) : RN_PI_OVER_TWO;
      if (normal_angle > max_pair_normal_angle) return 0;
      RNScalar normal_angle_factor = 1.0 - normal_angle / max_pair_normal_angle;
      affinity *= normal_angle_factor;
    }
    else if ((primitive.primitive_type == LINE_PRIMITIVE_TYPE) && (cluster->primitive.primitive_type == PLANE_PRIMITIVE_TYPE)) {
      RNScalar dot = fabs(primitive.line.Vector().Dot(cluster->primitive.plane.Normal()));
      RNAngle normal_angle = (dot < 1) ? RN_PI_OVER_TWO - acos(dot) : RN_PI_OVER_TWO;
      if (normal_angle > max_pair_normal_angle) return 0;
      RNScalar normal_angle_factor = 1.0 - normal_angle / max_pair_normal_angle;
      affinity *= normal_angle_factor;
    }
    else if ((primitive.primitive_type == PLANE_PRIMITIVE_TYPE) && (cluster->primitive.primitive_type == PLANE_PRIMITIVE_TYPE)) {
      RNScalar dot = fabs(primitive.plane.Normal().Dot(cluster->primitive.plane.Normal()));
      RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
      if (normal_angle > max_pair_normal_angle) return 0;
      RNScalar normal_angle_factor = 1.0 - normal_angle / max_pair_normal_angle;
      affinity *= normal_angle_factor;
    }
  }

#if 0
  // Compute imbalance
  double min_imbalance = 0.01;
  double max_imbalance = 0;
  if ((min_imbalance > 0) || (max_imbalance > 0)) {
    int npoints1 = (points.NEntries() < cluster->points.NEntries()) ? points.NEntries() : cluster->points.NEntries();
    int npoints2 = (points.NEntries() > cluster->points.NEntries()) ? points.NEntries() : cluster->points.NEntries();
    if (npoints2 == 0) return 0;
    RNScalar imbalance = (double) npoints1 / (double) npoints2;
    if ((max_imbalance > 0) && (imbalance > max_imbalance)) return 0;
    if (imbalance < min_imbalance) imbalance = min_imbalance;
    affinity *= imbalance;
  }
#endif
  
  // Return affinity
  return affinity;
}



static int
CompareClusters(const void *data1, const void *data2)
{
  Cluster *cluster1 = *((Cluster **) data1);
  Cluster *cluster2 = *((Cluster **) data2);
  if (cluster2->total_affinity > cluster1->total_affinity) return 1;
  else if (cluster1->total_affinity > cluster2->total_affinity) return -1;
  else return 0;
}



////////////////////////////////////////////////////////////////////////
// Pair member functions
////////////////////////////////////////////////////////////////////////

Pair::
Pair(Cluster *cluster1, Cluster *cluster2, RNScalar affinity)
  : affinity(affinity),
    heapentry(NULL)
{
  // Insert pair into clusters
  if (cluster1 && cluster2) {
    // Remember clusters
    clusters[0] = cluster1;
    clusters[1] = cluster2;

    // Remember position of pair in clusters
    cluster_index[0] = cluster1->pairs.NEntries();
    cluster_index[1] = cluster2->pairs.NEntries();

    // Update clusters
    cluster1->pairs.Insert(this);
    cluster2->pairs.Insert(this);
  }
  else {
    // Initialize clusters
    clusters[0] = NULL;
    clusters[1] = NULL;

    // Initialize cluster index
    cluster_index[0] = -1;
    cluster_index[1] = -1;
  }
}



Pair::
~Pair(void)
{
  // Remove this pair from first cluster
  if (clusters[0]) {
    assert(cluster_index[0] >= 0);
    RNArrayEntry *entry = clusters[0]->pairs.KthEntry(cluster_index[0]);
    Pair *tail = clusters[0]->pairs.Tail();
    if (tail->clusters[0] == clusters[0]) tail->cluster_index[0] = cluster_index[0];
    else if (tail->clusters[1] == clusters[0]) tail->cluster_index[1] = cluster_index[0];
    clusters[0]->pairs.EntryContents(entry) = tail;
    clusters[0]->pairs.RemoveTail();
  }

  // Remove this pair from second cluster
  if (clusters[1]) {
    assert(cluster_index[1] >= 0);
    RNArrayEntry *entry = clusters[1]->pairs.KthEntry(cluster_index[1]);
    Pair *tail = clusters[1]->pairs.Tail();
    if (tail->clusters[0] == clusters[1]) tail->cluster_index[0] = cluster_index[1];
    else if (tail->clusters[1] == clusters[1]) tail->cluster_index[1] = cluster_index[1];
    clusters[1]->pairs.EntryContents(entry) = tail;
    clusters[1]->pairs.RemoveTail();
  }
}



static Pair *
FindPair(Cluster *cluster1, Cluster *cluster2) 
{
  // Swap clusters so that cluster1 has fewer pairs
  if (cluster1->pairs.NEntries() > cluster2->pairs.NEntries()) {
    Cluster *swap = cluster1; 
    cluster1 = cluster2; 
    cluster2 = swap;
  }

  // Search for pair
  for (int i = 0; i < cluster1->pairs.NEntries(); i++) {
    Pair *pair = cluster1->pairs.Kth(i);
    if (pair->clusters[0] == cluster2) return pair;
    if (pair->clusters[1] == cluster2) return pair;
  }

  // Pair not found
  return NULL;
}



////////////////////////////////////////////////////////////////////////
// Segmentation functions
////////////////////////////////////////////////////////////////////////

Segmentation::
Segmentation(void)
  : points(),
    kdtree(NULL),
    clusters(),
    point_buffer(NULL)
{
}



Segmentation::
~Segmentation(void)
{
  // Delete clusters
  for (int i = 0; i < clusters.NEntries(); i++) delete clusters[i];
  
  // Delete kdtree
  if (kdtree) delete kdtree;

  // Delete points
  if (point_buffer) delete [] point_buffer;
  else { for (int i = 0; i < points.NEntries(); i++) delete points[i]; }
}



RNScalar Segmentation::
Affinity(void) const
{
  RNScalar sum = 0;
  for (int i = 0; i < clusters.NEntries(); i++) {
    sum += clusters[i]->total_affinity;
  }
  return sum;
}



int Segmentation::
NUnclusteredPoints(void) const
{
  // Count unclustered points
  int count = 0;
  for (int i = 0; i < points.NEntries(); i++) {
    Point *point = points.Kth(i);
    if (!point->cluster) count++;
  }

  // Return number of unclustered points
  return count;
}



int Segmentation::
CreatePoints(const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image, 
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image,
  const R2Grid& boundary_image, const R2Image& color_image)
{
  // Allocate points 
  point_buffer = new Point [ depth_image.NEntries() ];
  if (!point_buffer) {
    fprintf(stderr, "Unable to allocate points\n");
    return 0;
  }

  // Fill points
  for (int ix = 0; ix < depth_image.XResolution(); ix++) {
    for (int iy = 0; iy < depth_image.YResolution(); iy++) {
      int i;
      depth_image.IndicesToIndex(ix, iy, i);
      Point *point = &point_buffer[i];

      // Check depth
      RNScalar depth = depth_image.GridValue(i);
      if (RNIsNegativeOrZero(depth)) continue;
      if ((max_depth > 0) && (depth > max_depth)) continue;
      point->depth = depth;

      // Get position
      RNScalar px = px_image.GridValue(i);
      RNScalar py = py_image.GridValue(i);
      RNScalar pz = pz_image.GridValue(i);
      point->position.Reset(px, py, pz);

      // Get normal
      RNScalar nx = nx_image.GridValue(i);
      RNScalar ny = ny_image.GridValue(i);
      RNScalar nz = nz_image.GridValue(i);
      point->normal.Reset(nx, ny, nz);

      // Get radius
      RNScalar radius = radius_image.GridValue(i);
      point->radius = radius;

      // Get color
      point->color = color_image.PixelRGB(ix, iy);
    
      // Get flags
      point->boundary = (unsigned int) (boundary_image.GridValue(i) + 0.5);

      // Set grid index
      point->grid_index = i;

      // Insert point
      points.Insert(point);
    }
  }

  // Create kdtree of points
  Point tmp; int position_offset = (unsigned char *) &(tmp.position) - (unsigned char *) &tmp;
  kdtree = new R3Kdtree<Point *>(points, position_offset);
  if (!kdtree) {
    fprintf(stderr, "Unable to create kdtree\n");
    return 0;
  }
  
  // Create arrays of neighbor points
  for (int i = 0; i < points.NEntries(); i++) {
    Point *point = points.Kth(i);
    int ix, iy, neighbor_index;
    depth_image.IndexToIndices(point->grid_index, ix, iy);
    for (int s = -1; s <= 1; s++) {
      if ((ix+s < 0) || (ix+s >= depth_image.XResolution())) continue;
      for (int t = -1; t <= 1; t++) {
        if ((s == 0) && (t == 0)) continue;
        if ((iy+t < 0) || (iy+t >= depth_image.YResolution())) continue;
        depth_image.IndicesToIndex(ix+s, iy+t, neighbor_index);
        Point *neighbor = &point_buffer[neighbor_index];
        if ((point->boundary & FEATURE_IS_ON_SHADOW_BOUNDARY) && (neighbor->boundary & FEATURE_IS_ON_SILHOUETTE_BOUNDARY)) continue;
        if ((point->boundary & FEATURE_IS_ON_SILHOUETTE_BOUNDARY) && (neighbor->boundary & FEATURE_IS_ON_SHADOW_BOUNDARY)) continue;
        point->neighbors.Insert(neighbor);
      }
    }
  }

  // Return success
  return 1;
}
               


int Segmentation::
CreateSingletonClusters(int primitive_type)
{
  // Create cluster for every point
  for (int i = 0; i < points.NEntries(); i++) {
    Point *point = points.Kth(i);

    // Create primitive
    Primitive primitive(primitive_type);
    primitive.Update(point);
    
    // Create cluster
    Cluster *cluster = new Cluster(point, primitive);

    // Insert point
    cluster->InsertPoint(point, 1.0);

    // Insert cluster
    cluster->segmentation_index = clusters.NEntries();    
    clusters.Insert(cluster);
  }

  // Return success
  return 1;
}



int Segmentation::
CreateRansacClusters(int primitive_type)
{
  // Check number of ransac iterations
  if (max_ransac_iterations == 0) return 1;
  
  // Determine how many seed points to skip each iteration
  int skip = 1;
  if ((max_clusters > 0) && (points.NEntries()/(4*max_clusters) > skip))
    skip = points.NEntries()/(4*max_clusters);
  if ((min_cluster_points > 0) & ((min_cluster_points/4) > skip))
    skip = min_cluster_points/4;
  if ((min_clusters > 0) && (skip > points.NEntries()/min_clusters))
    skip = points.NEntries()/min_clusters;

  // Search seed points
  int seed_index = 0;
  while (seed_index < points.NEntries()) {
    // Find next seed point
    Point *seed_point = NULL;
    while ((seed_index < points.NEntries()) && !seed_point) {
      Point *point = points.Kth(seed_index);
      if (!point->cluster || (point->cluster_affinity < 0.1)) seed_point = point;
      seed_index += skip; 
    }

    // Check seed point
    if (!seed_point) break;

    // Create cluster
    Primitive primitive(primitive_type);
    primitive.Update(seed_point);
    Cluster *cluster = new Cluster(seed_point, primitive);
    if (!cluster->UpdatePoints(kdtree)) { delete cluster; continue; }

    // Iteratively update everything
    RNBoolean error = FALSE;
    for (int iter = 0; iter < max_ransac_iterations; iter++) {
      if (!cluster->UpdatePrimitive()) { error = TRUE; break; }
      if (!cluster->UpdatePoints(kdtree)) { error = TRUE; break; }
    }

    // Check for error
    if (error) { 
      delete cluster; 
      continue; 
    }

    // Insert cluster
    cluster->segmentation_index = clusters.NEntries();    
    clusters.Insert(cluster);
  } 

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Clustering manipulation functions
////////////////////////////////////////////////////////////////////////

int Segmentation::
RefineClusters(void)
{
  // Iteratively update everything
  int max_iterations = 1;
  for (int iter = 0; iter < max_iterations; iter++) {
    // Copy list of clusters
    RNArray<Cluster *> tmp = clusters;
    tmp.Sort(CompareClusters);

    // Will rebuild list of clusters
    clusters.Empty();

    // Refine each cluster
    RNBoolean converged = TRUE;
    for (int i = 0; i < tmp.NEntries(); i++) {
      Cluster *cluster = tmp.Kth(i);
      int prev_npoints = cluster->points.NEntries();

      // Refine cluster
      RNBoolean error = FALSE;
      if (!error && !cluster->UpdatePrimitive()) error = TRUE;
      if (!error && !cluster->UpdatePoints(kdtree)) error = TRUE; 

      // Insert cluster
      cluster->segmentation_index = clusters.NEntries();    
      if (!error) clusters.Insert(cluster);
      else delete cluster;

      // Check for convergence
      if (error || (prev_npoints != cluster->points.NEntries())) {
        converged = FALSE;
      }
    }

    // Check if converged
    if (converged) break;
  }

  // Return success
  return 1;
}



int Segmentation::
DeleteClusters(void)
{
  // Sort clusters
  clusters.Sort(CompareClusters);

  // Separate viable from nonviable ones
  RNArray<Cluster *> viable_clusters;
  RNArray<Cluster *> nonviable_clusters;
  for (int i = 0; i < clusters.NEntries(); i++) {
    Cluster *cluster = clusters.Kth(i);

    // Check min_clusters
    if ((min_clusters <= 0) || (i >= min_clusters)) {
      // Check cluster points
      if (min_cluster_points > 0) {
        if (cluster->points.NEntries() < min_cluster_points) {
          nonviable_clusters.Insert(cluster);
          continue;
        }
      }

      // Check cluster coverage
      if (min_cluster_coverage > 0) {
        if (cluster->Coverage() < min_cluster_coverage) {
          nonviable_clusters.Insert(cluster);
          continue;
        }
      }

      // Check max_clusters
      if (max_clusters > 0) {
        if (viable_clusters.NEntries() > max_clusters) {
          nonviable_clusters.Insert(cluster);
          continue;
        }
      }
    }
    
    // Cluster is viable
    cluster->segmentation_index = viable_clusters.NEntries();    
    viable_clusters.Insert(cluster);
  }

  // Delete nonviable clusters
  for (int i = 0; i < nonviable_clusters.NEntries(); i++) {
    Cluster *cluster = nonviable_clusters.Kth(i);
    delete cluster;
  }

  // Replace clusters with viable ones
  clusters = viable_clusters;

  // Return success
  return 1;
}



int Segmentation::
MergeClusters(void)
{
  // Initialize statistics
  int merge_count = 0;
  int push_count = 0;

  //////////

  // Create pairs between clusters with nearby points
  RNArray<Pair *> pairs;
  for (int i = 0; i < clusters.NEntries(); i++) {
    Cluster *cluster0 = clusters.Kth(i);
    if (cluster0->points.IsEmpty()) continue;
    Point *point0 = cluster0->points.Head();

    // Create pairs
    for (int j = 0; j < point0->neighbors.NEntries(); j++) {
      Point *point1 = point0->neighbors.Kth(j);
      if (point0 == point1) continue;
      Cluster *cluster1 = point1->cluster;
      if (!cluster1) continue;
      if (cluster0 == cluster1) continue;

      // Check if within max neighbor distance
      if (max_neighbor_distance_factor > 0) {
        RNScalar radius = (point0->radius > point1->radius) ? point1->radius : point0->radius;
        RNScalar max_neighbor_distance = (radius > 0) ? max_neighbor_distance_factor * radius : 0.25;
        RNScalar dd = R3SquaredDistance(point0->position, point1->position);
        if (dd > max_neighbor_distance * max_neighbor_distance) continue;
      }

      // Check if already have pair
      if (FindPair(cluster0, cluster1)) continue;

      // Compute affinity
      RNScalar affinity = cluster0->Affinity(cluster1);
      if (affinity < min_pair_affinity) continue;

      // Create pair
      Pair *pair = new Pair(cluster0, cluster1, affinity);
      if (!pair) continue;

      // Insert pair
      pairs.Insert(pair);
    }
  }

  // Check if there are any pairs
  if (pairs.IsEmpty()) return 1;

  //////////

  // Initialize heap
  Pair tmp;
  RNHeap<Pair *> heap(&tmp, &tmp.affinity, &tmp.heapentry, FALSE);
  for (int i = 0; i < pairs.NEntries(); i++) {
    Pair *pair = pairs.Kth(i);
    heap.Push(pair);
  }

  // Merge clusters hierarchically
  while (!heap.IsEmpty()) {
    // Get pair
    Pair *pair = heap.Pop();

    // Check if we are done
    if (pair->affinity < min_pair_affinity) break;

    // Get clusters
    Cluster *cluster0 = pair->clusters[0];
    Cluster *cluster1 = pair->clusters[1];

    // Check if either cluster has already been merged
    if (cluster0->parent || cluster1->parent) {
      // Find ancestors
      Cluster *ancestor0 = cluster0;
      Cluster *ancestor1 = cluster1;
      while (ancestor0->parent) ancestor0 = ancestor0->parent;
      while (ancestor1->parent) ancestor1 = ancestor1->parent;
      if (ancestor0 != ancestor1) {
        if (!FindPair(ancestor0, ancestor1)) {
          RNScalar affinity = ancestor0->Affinity(ancestor1);
          if (affinity > min_pair_affinity) {
            // Create a pair between the ancestors
            Pair *pair = new Pair(ancestor0, ancestor1, affinity);
            heap.Push(pair);
            push_count++;
          }
        }
      }
    }
    else {
      if (0 && print_progress) {
        static unsigned long count = 0;
        if ((count++ % 1000) == 0) {
          printf("        %15.12f : %9d %9d : %15d %15d %15d\n", pair->affinity, 
                 cluster0->points.NEntries(), cluster1->points.NEntries(), 
                 heap.NEntries(), merge_count, push_count);
        }
      }

#if 0
      // Create merged cluster
      Cluster *cluster = new Cluster(cluster0, cluster1);
      clusters.Insert(cluster);
      merge_count++;
#else
      // Merge smaller cluster into bigger one
      Cluster *parent = (cluster0->points.NEntries() > cluster1->points.NEntries()) ? cluster0 : cluster1;
      Cluster *child = (cluster0->points.NEntries() > cluster1->points.NEntries()) ? cluster1 : cluster0;
      parent->InsertChild(child);
      merge_count++;
#endif
    }

    // Delete pair
    delete pair;
  }

  // Remove merged clusters
  RNArray<Cluster *> merged_clusters;
  RNArray<Cluster *> all_clusters = clusters;
  clusters.Empty();
  for (int i = 0; i < all_clusters.NEntries(); i++) {
    Cluster *cluster = all_clusters.Kth(i);
    cluster->segmentation_index = clusters.NEntries();    
    if (!cluster->parent) { clusters.Insert(cluster); continue; }
    cluster->parent->RemoveChild(cluster);
    merged_clusters.Insert(cluster);
  }

  // Delete merged clusters
  for (int i = 0; i < merged_clusters.NEntries(); i++) {
    Cluster *cluster = merged_clusters.Kth(i);
    delete cluster;
  }

  // Return success
  return 1;
}



int Segmentation::
SplitClusters(void)
{
#if 0
  // Check min clusters spacing
  if (min_cluster_spacing <= 0) return 1;

  // Split connected components
  RNArray<Cluster *> tmp = clusters;
  clusters.Empty();
  for (int i = 0; i < tmp.NEntries(); i++) {
    Cluster *cluster = tmp.Kth(i);

    // Check cluster
    if (cluster->points.NEntries() < min_cluster_points) continue;

    // Rasterize points into planar grid
    R3PlanarGrid grid(cluster->primitive.plane, cluster->primitive.bbox, min_cluster_spacing);
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      Point *point = cluster->points.Kth(j);
      R3Point position = point->position;
      grid.RasterizeWorldPoint(position.X(), position.Y(), position.Z(), 1.0);
    }

    // Compute connected components
    int max_components = grid.NEntries();
    int *components = new int [ max_components ];
    int ncomponents = grid.ConnectedComponents(RN_EPSILON, max_components, NULL, NULL, components);

    // Check connected components
    if (ncomponents == 1) {
      // One connected component - simply insert cluster
      cluster->segmentation_index = clusters.NEntries();    
      clusters.Insert(cluster);
    }
    else {
      // Create cluster for each connnected component
      for (int j = 0; j < ncomponents; j++) {
        // Make array of points in component
        RNArray<Point *> component_points;
        for (int k = 0; k < cluster->points.NEntries(); k++) {
          Point *point = cluster->points.Kth(k);
          R3Point world_position = point->position;
          R2Point grid_position = grid.GridPosition(world_position);
          int ix = (int) (grid_position.X() + 0.5);
          int iy = (int) (grid_position.Y() + 0.5);
          int index; grid.Grid().IndicesToIndex(ix, iy, index);
          if (components[index] != j) continue;
          component_points.Insert(point);
        }

        // Check number of points
        if (component_points.NEntries() > min_cluster_points) {

          // Find centroid
          R3Point centroid = R3zero_point;
          for (int k = 0; k < component_points.NEntries(); k++) {
            Point *point = component_points.Kth(k);
            R3Point world_position = point->position;
            centroid += world_position;
          }
          centroid /= component_points.NEntries();

          // Find seed point
          Point *seed_point = NULL;
          RNLength min_dd = FLT_MAX;
          for (int k = 0; k < component_points.NEntries(); k++) {
            Point *point = component_points.Kth(k);
            R3Point world_position = point->position;
            RNLength dd = R3SquaredDistance(centroid, world_position);
            if (dd < min_dd) { seed_point = point; min_dd = dd; }
          }

          // Check seed point
          if (seed_point) {
            // Create cluster
            Cluster *c = new Cluster(seed_point, PLANE_PRIMITIVE_TYPE);
            c->possible_affinity = cluster->possible_affinity;

            // Insert points into cluster
            for (int k = 0; k < component_points.NEntries(); k++) {
              Point *point = component_points.Kth(k);
              c->InsertPoint(point);
            }

            // Update primitive
            c->UpdatePrimitive();

            // Update planar grid
            // c->UpdatePlanarGrid();

            // Insert cluster
            c->segmentation_index = clusters.NEntries();    
            clusters.Insert(c);
          }
        }
      }

      // Delete the original cluster
      delete cluster;
    }

    // Delete components
    delete [] components;
  }
#endif

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Top-level segmentation functions
////////////////////////////////////////////////////////////////////////

int Segmentation::
CreateClusters(int primitive_type)
{
  // Print debug message
  RNTime step_time;
  step_time.Read();
  if (print_progress) {
    printf("      SA %.3f %d\n", step_time.Elapsed(), points.NEntries());
    step_time.Read();
  }

  // Create clusters
  if (initialize_hierarchically) {
    if (!CreateSingletonClusters(primitive_type)) return 0;
    if (!MergeClusters()) return 0;
  }
  else {
    if (!CreateRansacClusters(primitive_type)) return 0;
  }

  // Check clusters
  if (clusters.IsEmpty()) return 0;

  // Print debug message
  if (print_progress) {
    printf("      SB %.3f %d %d %g\n", step_time.Elapsed(), clusters.NEntries(), NUnclusteredPoints(), Affinity());
    step_time.Read();
  }

  // Iteratively update clusters
  for (int i = 0; i < max_refinement_iterations; i++) {
    // Refine clusters
    if (!RefineClusters()) return 0;

    // Print debug message
    if (print_progress) {
      printf("      SC %d : %.3f %d %d %g\n", i, step_time.Elapsed(), clusters.NEntries(), NUnclusteredPoints(), Affinity());
      step_time.Read();
    }

    // Create clusters
    if (!CreateRansacClusters(primitive_type)) return 0;

    // Print debug message
    if (print_progress) {
      printf("      SD %d : %.3f %d %d %g\n", i, step_time.Elapsed(), clusters.NEntries(), NUnclusteredPoints(), Affinity());
      step_time.Read();
    }

    // Merge clusters
    if (!MergeClusters()) return 0;

    // Print debug message
    if (print_progress) {
      printf("      SE %d : %.3f %d %d %g\n", i, step_time.Elapsed(), clusters.NEntries(), NUnclusteredPoints(), Affinity());
      step_time.Read();
    }

    // Delete clusters
    if (!DeleteClusters()) return 0;

    // Print debug message
    if (print_progress) {
      printf("      SF %d : %.3f %d %d %g\n", i, step_time.Elapsed(), clusters.NEntries(), NUnclusteredPoints(), Affinity());
      step_time.Read();
    }
  }

  // Split clusters
  // if (!SplitClusters()) return 0;

  // Print debug message
  if (print_progress) {
    printf("      SG %.3f %d %d %g\n", step_time.Elapsed(), clusters.NEntries(), NUnclusteredPoints(), Affinity());
    step_time.Read();
  }

  // Delete clusters
  if (!DeleteClusters()) return 0;

  // Print debug message
  if (print_progress) {
    printf("      SH %.3f %d %d %g\n", step_time.Elapsed(), clusters.NEntries(), NUnclusteredPoints(), Affinity());
    step_time.Read();
  }

  // Sort clusters
  clusters.Sort(CompareClusters);

  // Return success
  return 1;
}




int Segmentation::
ReadClusterImage(const char *filename)
{
  // Read file
  R2Grid segmentation_grid;
  if (!segmentation_grid.ReadFile(filename)) return 0;

  // Process segmentation grid
  segmentation_grid.Substitute(0, R2_GRID_UNKNOWN_VALUE);
  segmentation_grid.Subtract(1.0);

  // Create clusters
  for (int i = 0; i < segmentation_grid.NEntries(); i++) {
    Point *point = &point_buffer[i];
    if (point->depth < 0) continue;
    RNScalar grid_value = segmentation_grid.GridValue(i);
    if (grid_value == R2_GRID_UNKNOWN_VALUE) continue;
    int cluster_index = (int) (grid_value + 0.5);
    while (clusters.NEntries() <= cluster_index) {
      // Create cluster
      Cluster *cluster = new Cluster(NULL, PLANE_PRIMITIVE_TYPE);
      cluster->segmentation_index = clusters.NEntries();    
      clusters.Insert(cluster);
    }
    Cluster *cluster = clusters.Kth(cluster_index);
    if (!cluster->seed_point) cluster->seed_point = point;
    cluster->InsertPoint(point);
  }

  // Update cluster primitives
  for (int i = 0; i < clusters.NEntries(); i++) {
    Cluster *cluster = clusters.Kth(i);
    cluster->UpdatePrimitive();
  }
  
  // Sort clusters
  clusters.Sort(CompareClusters);

  // Return success
  return 1;
}



int Segmentation::
WriteClusterImage(int xres, int yres, const char *filename) const
{
  // Fill image
  R2Grid image(xres, yres);
  for (int i = 0; i < clusters.NEntries(); i++) {
    Cluster *cluster = clusters.Kth(i);
    for (int j = 0; j < cluster->points.NEntries(); j++) {
      Point *point = cluster->points.Kth(j);
      if (point->grid_index < 0) continue;
      image.SetGridValue(point->grid_index, i+1);
    }
  }

  // Write image
  if (!image.WriteFile(filename)) return 0;
  
  // Return success
  return 1;   
}



////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////

static int
FindPairIndex(RNArray<Pair *>& pairs, Cluster *cluster0, Cluster *cluster1)
{
  // Search for primitive marker
  for (int i = 0; i < pairs.NEntries(); i++) {
    Pair *pair = pairs.Kth(i);
    if ((cluster0 == pair->clusters[0]) && (cluster1 == pair->clusters[1])) return i;
    else if ((cluster0 == pair->clusters[1]) && (cluster1 == pair->clusters[0])) return i;
  }

  // Not found
  return -1;
}



static Segmentation *
CreateSegmentation(const R2Grid& px_image, const R2Grid& py_image, const R2Grid& pz_image,
  const R2Grid& nx_image, const R2Grid& ny_image, const R2Grid& nz_image,
  const R2Grid& depth_image, const R2Grid& radius_image, 
  const R2Grid& boundary_image, const R2Image& color_image,
  const R3Point& viewpoint, const R3Vector& towards, const R3Vector& up,
  const char *segmentation_filename = NULL)
{
  // Adjust segmentation parameters ???
  min_cluster_points = 10 * depth_image.NEntries() / (640 * 480);
  
  // Allocate segmentation
  Segmentation *segmentation = new Segmentation();
  if (!segmentation) {
    fprintf(stderr, "Unable to allocate segmentation.\n");
    return NULL;
  }

  // Create points
  if (!segmentation->CreatePoints(px_image, py_image, pz_image, nx_image, ny_image, nz_image,
    depth_image, radius_image, boundary_image, color_image)) {
    fprintf(stderr, "Unable to create points for segmentation.\n");
    delete segmentation;
    return 0;
  }

  // Check points
  if (segmentation->points.NEntries() == 0) {
    delete segmentation;
    return 0;
  }

  // Check if filename given
  if (segmentation_filename && strcmp(segmentation_filename, "-")) {
    // Read clusters
    if (!segmentation->ReadClusterImage(segmentation_filename)) {
      delete segmentation;
      return 0;
    }
  }
  else {
    // Create clusters
    if (!segmentation->CreateClusters(PLANE_PRIMITIVE_TYPE)) {
      delete segmentation;
      return 0;
    }
  }

  // Return segmentation
  return segmentation;
}



