// FineToCoarseRegistration - fetregister - planar proxy clustering
//
// This file contains implementation of agglomerative clustering algorithm to 
// create structural model, namely set of base proxies P_b and cluster proxies 
// P_c.

// This code has been mostly written by Thomas Funkhouser.
// Minor modifications by Maciej Halber


////////////////////////////////////////////////////////////////////////////////
// License:
//
// Copyright (c) 2017 Maciej Halber, Thomas Funkhouser
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

#ifndef STRUCTURE_H_
#define STRUCTURE_H_


////////////////////////////////////////////////////////////////////////////////
// Type definitions
////////////////////////////////////////////////////////////////////////////////

struct Proxy {
public:
  Proxy(struct Structure *structure);
  Proxy(struct Structure *structure, const RNArray<FETFeature *>& features, int shape_index);
  Proxy(struct Structure *structure, Proxy *child);
  Proxy(struct Structure *structure, Proxy *child1, Proxy *child2);
  ~Proxy(void);

  void InsertChild(Proxy *proxy);
  void RemoveChild(Proxy *proxy);

public:
  // Structure membership
  struct Structure *structure;
  int structure_index;

  // Proxy hierarchy
  Proxy *parent;
  RNArray<Proxy *> children;

  // Fet stuff
  FETShape *shape;
  FETFeature *feature;

  // Hierarchical clustering stuff
  RNArray<struct Pair *> pairs;
  Proxy *cluster_representative;

  // "Segment" stuff
  int min_shape_index;
  int max_shape_index;
  int level;

  // Properties
  int n_inliers;
  double mean;
  double sdev;
};

struct Pair {
public:
  Pair(Proxy *proxy1 = NULL, Proxy *proxy2 = NULL, RNScalar affinity = 0);
  ~Pair(void);
public:
  Proxy *proxies[2];
  int proxy_index[2];
  RNScalar affinity;
  Pair **heapentry;
};

struct Structure {
public:
  Structure(FETReconstruction *reconstruction = NULL);
  ~Structure(void);

  void InsertProxy(Proxy *proxy);
  void RemoveProxy(Proxy *proxy);
  void CleanProxies(void);

  Proxy * GetProxy( int k ) const;
  int NProxies() const;

  int CreateBaseLevelProxies( double max_depth );
  int CreateNextLevelProxies( int level, 
                              double *parametrization, double segment_length);

public:
  FETReconstruction *reconstruction;
  RNArray<Proxy *> proxies;
  R3Box bbox;
};


////////////////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////////////

static int min_features_per_proxy = 100;
static int min_features_per_cluster = 2000;
static double min_proxy_normal_mean = 0.3;
static double max_pair_centroid_distance = 256;
static double max_pair_plane_distance = 0.3;
static double max_pair_normal_angle = RN_PI / 8.0;
static double min_pair_affinity = 1.0E-6;

////////////////////////////////////////////////////////////////////////////////
// Proxy member functions
////////////////////////////////////////////////////////////////////////////////

Proxy::
Proxy(Structure *structure)
  : structure(NULL),
    structure_index(-1),
    parent(NULL),
    children(),
    shape(NULL),
    feature(NULL),
    pairs(),
    cluster_representative(NULL),
    min_shape_index(-1),
    max_shape_index(-1),
    level(-1),
    n_inliers(0)
{
  // Insert proxy into structure
  structure->InsertProxy(this);

  // Create shape
  this->shape = new FETShape(structure->reconstruction);

  // Create feature
  this->feature = new FETFeature(structure->reconstruction, PLANE_FEATURE_SHAPE);
  feature->generator_type = STRUCTURE_FEATURE_TYPE;
  this->shape->InsertFeature(this->feature);
}



Proxy::
Proxy(Structure *structure, const RNArray<FETFeature *>& features, int shape_index)
  : structure(NULL),
    structure_index(-1),
    parent(NULL),
    children(),
    shape(NULL),
    feature(NULL),
    pairs(),
    cluster_representative(NULL),
    min_shape_index(-1),
    max_shape_index(-1),
    level(-1),
    n_inliers(0)
{
  // Insert proxy into structure
  structure->InsertProxy(this);

  // Create shape
  this->shape = new FETShape(structure->reconstruction);

  // Create feature
  this->feature = new FETFeature(structure->reconstruction, PLANE_FEATURE_SHAPE);
  feature->generator_type = STRUCTURE_FEATURE_TYPE;
  this->shape->InsertFeature(this->feature);

  // Check features
  if (features.IsEmpty()) return;

  // Compute feature properties
  RNScalar salience = 0;
  R3Point position = R3zero_point;
  R3Vector normal = R3zero_vector;
  for (int i = 0; i < features.NEntries(); i++) {
    FETFeature *f = features.Kth(i);
    salience += f->Salience();
    position += f->Position(TRUE);
    normal += f->Normal(TRUE);
  }

  // Set feature properties
  if (salience > 0) {
    // Compute feature position
    position /= features.NEntries();

    // Compute feature radius
    RNScalar radius = 0;
    for (int i = 0; i < features.NEntries(); i++) {
      FETFeature *f = features.Kth(i);
      RNScalar d = R3Distance(f->Position(TRUE), position);
      radius += f->Salience() * d;
    }

    // Set properties
    normal /= features.NEntries();
    normal.Normalize();
    radius /= salience;
    feature->SetPosition(position);
    feature->SetNormal(normal);
    feature->SetSalience(salience);
    feature->SetRadius(radius);
    feature->SetColor(RNred_rgb);
  }

  // Create correspondences
  for (int i = 0; i < features.NEntries(); i++) {
    FETFeature *f = features.Kth(i);
    new FETCorrespondence(structure->reconstruction, f, this->feature, 1.0, COINCIDENT_RELATIONSHIP);
  }

  // Set shape index range
  this->min_shape_index = shape_index;
  this->max_shape_index = shape_index;

  // Set level
  this->level = 0;
  this->n_inliers = features.NEntries();
}



Proxy::
Proxy(Structure *structure, Proxy *child)
  : structure(NULL),
    structure_index(-1),
    parent(NULL),
    children(),
    shape(NULL),
    feature(NULL),
    pairs(),
    cluster_representative(NULL),
    min_shape_index(-1),
    max_shape_index(-1),
    level(-1),
    n_inliers(0)
{
  // Insert proxy into structure
  structure->InsertProxy(this);

  // Create shape
  this->shape = new FETShape(structure->reconstruction);

  // Create feature
  this->feature = new FETFeature(*(child->feature));
  this->feature->SetPosition ( child->feature->Position( TRUE ) );
  this->feature->SetNormal ( child->feature->Normal( TRUE ) );
  feature->generator_type = STRUCTURE_FEATURE_TYPE;
  this->shape->InsertFeature(this->feature);

  // Set shape index range
  this->min_shape_index = child->min_shape_index;
  this->max_shape_index = child->max_shape_index;

  // Set level
  this->level = child->level + 1;
  this->n_inliers += this->n_inliers + child->n_inliers;

  // Insert child into this proxy
  InsertChild(child);
}



Proxy::
Proxy(Structure *structure, Proxy *merge1, Proxy *merge2)
  : structure(NULL),
    structure_index(-1),
    parent(NULL),
    children(),
    shape(NULL),
    feature(NULL),
    pairs(),
    cluster_representative(NULL),
    min_shape_index(-1),
    max_shape_index(-1),
    level(-1),
    n_inliers(0)
{
  // Insert proxy into structure
  structure->InsertProxy(this);

  // Create shape
  this->shape = new FETShape(structure->reconstruction);

  // Create feature
  this->feature = new FETFeature(structure->reconstruction, PLANE_FEATURE_SHAPE);
  feature->generator_type = STRUCTURE_FEATURE_TYPE;
  this->shape->InsertFeature(this->feature);

  // Set feature properties
  RNScalar s1 = merge1->feature->Salience();
  RNScalar s2 = merge2->feature->Salience();
  RNScalar salience = s1 + s2;
  RNScalar radius = R3Distance(merge1->feature->Position(TRUE), merge2->feature->Position(TRUE));
  radius += (merge1->feature->Radius() < merge2->feature->Radius()) ? merge1->feature->Radius() : merge2->feature->Radius();
  if (salience > 0) {
    R3Point position = (s1 * merge1->feature->Position(TRUE) + s2 * merge2->feature->Position(TRUE)) / salience;
    R3Vector normal = s1 * merge1->feature->Normal(TRUE) + s2 * merge2->feature->Normal(TRUE);
    normal.Normalize();
    this->feature->SetPosition(position);
    this->feature->SetNormal(normal);
    this->feature->SetSalience(salience);
    this->feature->SetRadius(radius);
    this->feature->SetColor(RNgreen_rgb);
  }

  // Move children from merge1
  RNArray<Proxy *> children1 = merge1->children;
  for (int i = 0; i < children1.NEntries(); i++) {
    Proxy *child = children1.Kth(i);
    merge1->RemoveChild(child);
    this->InsertChild(child);
  }

  // Move children from merge2
  RNArray<Proxy *> children2 = merge2->children;
  for (int i = 0; i < children2.NEntries(); i++) {
    Proxy *child = children2.Kth(i);
    merge2->RemoveChild(child);
    this->InsertChild(child);
  }

  // Set shape index range
  this->min_shape_index = (merge1->min_shape_index < merge2->min_shape_index) ? merge1->min_shape_index : merge2->min_shape_index;
  this->max_shape_index = (merge1->max_shape_index > merge2->max_shape_index) ? merge1->max_shape_index : merge2->max_shape_index;

  // Assign level
  this->level = merge1->level;
  this->n_inliers = merge1->n_inliers + merge2->n_inliers;

  // Mark merge1 and merge2 as "inactive"
  merge1->min_shape_index = -1;
  merge1->max_shape_index = -1;
  merge2->min_shape_index = -1;
  merge2->max_shape_index = -1;
}



Proxy::
~Proxy(void)
{
  // Remove from parent
  if (this->parent) this->parent->RemoveChild(this);

  // Remove children
  while ( this->children.NEntries() > 0 ) 
  {
    Proxy *child = children.Kth(children.NEntries() - 1);
    this->RemoveChild(child);
  }
  
  // Delete shape
  if (this->shape) {
    delete this->shape; 
  }

  // Remove from structure
  if (this->structure) this->structure->RemoveProxy(this);
}



void Proxy::
InsertChild(Proxy *child)
{
  // Just checking
  assert(structure == child->structure);

  child->parent = this;

  // Create correspondence
  new FETCorrespondence(structure->reconstruction, child->feature, this->feature, 1.0, COINCIDENT_RELATIONSHIP );

  // Add to array of children
  children.Insert(child);
}



void Proxy::
RemoveChild(Proxy *child)
{
  // Just checking
  assert(structure == child->structure);

  child->parent = NULL;

  // Find correspondence
  if ( feature )
  {
    FETCorrespondence *correspondence = NULL;
    for (int i = 0; i < feature->NCorrespondences(); i++) {
      FETCorrespondence *c = feature->Correspondence(i);
      if (c->Feature(0) == child->feature ) { correspondence = c; break; }
    }

    // Delete correspondence
    assert(correspondence);
    delete correspondence;
  }
  // Remove from array of children
  children.Remove(child);
}



////////////////////////////////////////////////////////////////////////////////
// Pair member functions
////////////////////////////////////////////////////////////////////////////////

Pair::
Pair(Proxy *proxy1, Proxy *proxy2, RNScalar affinity)
  : affinity(affinity),
    heapentry(NULL)
{
  // Insert pair into proxies
  if (proxy1 && proxy2) {
    // Remember proxies
    proxies[0] = proxy1;
    proxies[1] = proxy2;

    // Remember position of pair in proxies
    proxy_index[0] = proxy1->pairs.NEntries();
    proxy_index[1] = proxy2->pairs.NEntries();

    // Update proxies
    proxy1->pairs.Insert(this);
    proxy2->pairs.Insert(this);
  }
  else {
    // Initialize proxies
    proxies[0] = NULL;
    proxies[1] = NULL;

    // Initialize proxy index
    proxy_index[0] = -1;
    proxy_index[1] = -1;
  }
}



Pair::
~Pair(void)
{
  // Remove this pair from first proxy
  if (proxies[0]) {
    assert(proxy_index[0] >= 0);
    RNArrayEntry *entry = proxies[0]->pairs.KthEntry(proxy_index[0]);
    Pair *tail = proxies[0]->pairs.Tail();
    if (tail->proxies[0] == proxies[0]) tail->proxy_index[0] = proxy_index[0];
    else if (tail->proxies[1] == proxies[0]) tail->proxy_index[1] = proxy_index[0];
    proxies[0]->pairs.EntryContents(entry) = tail;
    proxies[0]->pairs.RemoveTail();
  }

  // Remove this pair from second proxy
  if (proxies[1]) {
    assert(proxy_index[1] >= 0);
    RNArrayEntry *entry = proxies[1]->pairs.KthEntry(proxy_index[1]);
    Pair *tail = proxies[1]->pairs.Tail();
    if (tail->proxies[0] == proxies[1]) tail->proxy_index[0] = proxy_index[1];
    else if (tail->proxies[1] == proxies[1]) tail->proxy_index[1] = proxy_index[1];
    proxies[1]->pairs.EntryContents(entry) = tail;
    proxies[1]->pairs.RemoveTail();
  }
}



static Pair *
FindPair(Proxy *proxy1, Proxy *proxy2)
{
  // Swap proxies so that proxy1 has fewer pairs
  if (proxy1->pairs.NEntries() > proxy2->pairs.NEntries()) {
    Proxy *swap = proxy1;
    proxy1 = proxy2;
    proxy2 = swap;
  }

  // Search for pair
  for (int i = 0; i < proxy1->pairs.NEntries(); i++) {
    Pair *pair = proxy1->pairs.Kth(i);
    if (pair->proxies[0] == proxy2) return pair;
    if (pair->proxies[1] == proxy2) return pair;
  }

  // Pair not found
  return NULL;
}



static RNScalar
PairAffinity( Proxy *proxy1, Proxy *proxy2, 
              double *parametrization, double segment_length )
{
  // Initialize affinity
  RNScalar affinity = 1;

  // Get features
  FETFeature *feature1 = proxy1->feature;
  FETFeature *feature2 = proxy2->feature;


  // Compute parametrization distance
  if( parametrization )
  {
    int max_lim = std::max( proxy1->max_shape_index, proxy2->max_shape_index );
    int min_lim = std::min( proxy1->min_shape_index, proxy2->min_shape_index );
    double dist = parametrization[max_lim] -
                  parametrization[min_lim];
    if( dist >= segment_length ) return 0;
  }

  // Compute centroid distance
  if (max_pair_centroid_distance > 0) {
    RNLength centroid_distance = R3Distance( feature1->Position(TRUE), 
                                             feature2->Position(TRUE));
    if (centroid_distance > max_pair_centroid_distance) return 0;
    RNScalar centroid_distance_factor = 1.0 - centroid_distance / max_pair_centroid_distance;
    affinity *= centroid_distance_factor;
  }

  // Compute point-plane distances
  if (max_pair_plane_distance > 0) {
    // Compute plane1 distance
    R3Point pos1 = feature1->Position( TRUE );
    R3Vector norm1 = feature1->Normal( TRUE );
    R3Plane plane1( pos1, norm1 );
    RNLength plane1_distance = R3Distance( plane1, feature2->Position(TRUE) );
    if (plane1_distance > max_pair_plane_distance) return 0;
    RNScalar plane1_distance_factor = 1.0 - plane1_distance / max_pair_plane_distance;
    affinity *= plane1_distance_factor;

    // Compute plane2 distance
    R3Point pos2 = feature2->Position( TRUE );
    R3Vector norm2 = feature2->Normal( TRUE );
    R3Plane plane2( pos2, norm2 );
    RNLength plane2_distance = R3Distance( plane2, feature1->Position(TRUE) );
    if (plane2_distance > max_pair_plane_distance) return 0;
    RNScalar plane2_distance_factor = 1.0 - plane2_distance / max_pair_plane_distance;
    affinity *= plane2_distance_factor;
  }

  // Compute normal angle
  if (max_pair_normal_angle > 0) {
    R3Vector n1 = feature1->Normal(TRUE);
    R3Vector n2 = feature2->Normal(TRUE);
    
    n1.Normalize();
    n2.Normalize();
    if( n1.Length() == 0.0 || n2.Length() == 0.0 ) return 0;
    RNScalar dot = n1.Dot(n2);
    RNAngle normal_angle = (dot < 1) ? acos(dot) : 0;
    if (normal_angle > max_pair_normal_angle) return 0;
    RNScalar normal_angle_factor = 1.0 - normal_angle / max_pair_normal_angle;
    affinity *= normal_angle_factor;
  }

  // Return affinity
  return affinity;
}



////////////////////////////////////////////////////////////////////////////////
// Structure member functions
////////////////////////////////////////////////////////////////////////////////

Structure::
Structure(FETReconstruction *reconstruction)
  : reconstruction(reconstruction),
    proxies()
{
}



Structure::
~Structure(void)
{
  // Delete proxies
  std::vector<Proxy*> to_be_removed;
  for ( int i = 0 ; i < this->NProxies() ; ++i )
  {
    Proxy * proxy = proxies[ i ];
    to_be_removed.push_back( proxy );
  }
  for ( int i = 0 ; i < to_be_removed.size() ; ++i )
  {
    delete to_be_removed[ i ];
  }
}


void Structure::
InsertProxy(Proxy *proxy)
{
  // Just checking
  assert(proxy->structure_index == -1);
  assert(proxy->structure == NULL);

  // Insert proxy
  proxy->structure = this;
  proxy->structure_index = proxies.NEntries();
  proxies.Insert(proxy);
}


//NOTE: Smth wrong here!
void Structure::
RemoveProxy(Proxy *proxy)
{
  // Just checking
  assert(proxy->structure_index >= 0);
  assert(proxy->structure_index < proxies.NEntries());
  assert(proxy->structure == this);

  // Remove proxy
  RNArrayEntry *entry = proxies.KthEntry(proxy->structure_index);
  Proxy *tail = proxies.Tail();
  tail->structure_index = proxy->structure_index;
  proxies.EntryContents(entry) = tail;
  proxies.RemoveTail();
  proxy->structure_index = -1;
  proxy->structure = NULL;
}



void Structure::
CleanProxies(void)
{
  // Separate viable from nonviable proxies
  RNArray<Proxy *> viable_proxies;
  RNArray<Proxy *> nonviable_proxies;
  for (int i = 0; i < proxies.NEntries(); i++) {
    Proxy *proxy = proxies.Kth(i);

    // Check if top level
    if (!proxy->parent) {
      // Check number of children
      if (proxy->children.NEntries() <= 1) {
        nonviable_proxies.Insert(proxy);
        continue;
      }

      if ( proxy->level == 0 )
      {
        continue;
      }

      if (proxy->n_inliers < min_features_per_cluster )
      {
        nonviable_proxies.Insert(proxy);
        continue;
      }

      // Check salience
      if (proxy->feature->Salience() <= RN_EPSILON) {
        nonviable_proxies.Insert(proxy);
        continue;
      }

      // Check shape index range
      if ((proxy->min_shape_index < 0) || (proxy->max_shape_index < 0)) {
        nonviable_proxies.Insert(proxy);
        continue;
      }

      // Check if has been merged
      if (proxy->cluster_representative) {
        nonviable_proxies.Insert(proxy);
        continue;
      }

    }

    // Proxy is viable
    proxy->structure_index = viable_proxies.NEntries();
    viable_proxies.Insert(proxy);
  }

  // Delete nonviable proxies
  for (int i = 0; i < nonviable_proxies.NEntries(); i++) {
    Proxy *proxy = nonviable_proxies.Kth(i);
    delete proxy;
  }

  // Replace proxies with viable ones
  proxies = viable_proxies;
}

int Structure::
NProxies( void ) const
{
  return this->proxies.NEntries();
}

Proxy * Structure::
GetProxy( int k ) const
{
  return this->proxies.Kth( k );
}

int Structure::
CreateBaseLevelProxies( double max_depth )
{
  // Create proxies for every shape
  RNArray<FETShape *> base_shapes = reconstruction->shapes;
  for (int i = 0; i < base_shapes.NEntries(); i++) {
    FETShape *shape = base_shapes.Kth(i);

    // Create arrays of features associated with all primitive markers
    RNArray<RNArray<FETFeature *> *> primitive_features;
    for (int j = 0; j < shape->NFeatures(); j++) {
      FETFeature *feature = shape->Feature(j);
      if (feature->ShapeType() != PLANE_FEATURE_SHAPE) continue;
      int primitive_marker = feature->primitive_marker;
      if (primitive_marker < 0) continue;
      if ( RNIsLess(feature->Position().Z(), -max_depth) ) continue; // Exclude far away
      while (primitive_features.NEntries() <= primitive_marker)
        primitive_features.Insert(new RNArray<FETFeature *>());
      primitive_features[primitive_marker]->Insert(feature);
    }

    // Create proxy for each primitive marker
    for (int j = 0; j < primitive_features.NEntries(); j++) {
      RNArray<FETFeature *> *features = primitive_features.Kth(j);

      // calculate mean normal
      R3Vector mean_normal = R3zero_vector;
      for ( int k = 0 ; k < features->NEntries(); k++ )
      {
        FETFeature * f = features->Kth(k);
        mean_normal += f->Normal( TRUE );
      }
      mean_normal /= features->NEntries();
      mean_normal.Normalize();

      // compute mean and sdev
      double mean = 0.0, sdev = 0.0;
      for ( int k = 0 ; k < features->NEntries(); k++ )
      {
        R3Vector normal = features->Kth(k)->Normal( TRUE );
        mean += acos( normal.Dot(mean_normal)); 
      }
      mean /= features->NEntries();
      for ( int k = 0 ; k < features->NEntries(); k++ )
      {
        R3Vector normal = features->Kth(k)->Normal( TRUE );
        double diff = mean - acos( normal.Dot(mean_normal));
        sdev += diff * diff;
      }
      sdev = sqrt( sdev / features->NEntries() );

      if (features->NEntries() < min_features_per_proxy) continue;
      if ( RNIsGreater( mean, min_proxy_normal_mean ) ) continue;

      Proxy *p = new Proxy(this, *features, i);
      p->mean = mean;
      p->sdev = sdev;
    }

    // Delete arrays of features associated with all primitive markers
    for (int j = 0; j < primitive_features.NEntries(); j++) {
      if (primitive_features[j]) delete primitive_features[j];
    }
  }
  // Return success
  return 1;
}

int Structure::
CreateNextLevelProxies(int level, double *parametrization, double segment_length)
{
  // Create proxies
  RNArray<Proxy *> level_proxies;
  for (int i = 0; i < proxies.NEntries(); i++) {
    Proxy *proxy = proxies.Kth(i);
    if (proxy->parent) continue;
    if (proxy->level != level-1) continue;
    if (proxy->min_shape_index < 0) continue;
    if (proxy->max_shape_index < 0) continue;
    if (proxy->feature->Salience() == 0) continue;
    Proxy *parent = new Proxy(this, proxy);
    assert(parent->level == level);
    level_proxies.Insert(parent);
  }

  // Create pairs between proxies
  RNArray<Pair *> pairs;
  for (int i = 0; i < level_proxies.NEntries(); i++) {
    Proxy *proxy0 = level_proxies.Kth(i);
    for (int j = i+1; j < level_proxies.NEntries(); j++) {
      Proxy *proxy1 = level_proxies.Kth(j);

      // Compute affinity
      RNScalar affinity = PairAffinity(proxy0, proxy1, parametrization, segment_length);
      if (affinity < min_pair_affinity) continue;

      // Create pair
      Pair *pair = new Pair(proxy0, proxy1, affinity);
      if (!pair) continue;

      // Insert pair
      pairs.Insert(pair);
    }
  }

  // Check if there are any pairs
  if (pairs.IsEmpty())
  {
    CleanProxies();
    return 1;
  }

  // Initialize heap
  Pair tmp;
  RNHeap<Pair *> heap(&tmp, &tmp.affinity, &tmp.heapentry, FALSE);
  for (int i = 0; i < pairs.NEntries(); i++) {
    Pair *pair = pairs.Kth(i);
    heap.Push(pair);
  }

  // Merge proxies hierarchically
  while (!heap.IsEmpty()) {
    // Get pair
    Pair *pair = heap.Pop();

    // Check if we are done
    if (pair->affinity < min_pair_affinity) break;

    // Get proxies
    Proxy *proxy0 = pair->proxies[0];
    Proxy *proxy1 = pair->proxies[1];

    // Check if either proxy has already been merged
    if (proxy0->cluster_representative || proxy1->cluster_representative) {
      // Find ancestors
      Proxy *ancestor0 = proxy0;
      Proxy *ancestor1 = proxy1;
      while (ancestor0->cluster_representative) ancestor0 = ancestor0->cluster_representative;
      while (ancestor1->cluster_representative) ancestor1 = ancestor1->cluster_representative;
      if (ancestor0 != ancestor1) {
        if (!FindPair(ancestor0, ancestor1)) {
          RNScalar affinity = PairAffinity( ancestor0, ancestor1, 
                                            parametrization, segment_length );
          if (affinity > min_pair_affinity) {
            // Create a pair between the ancestors
            Pair *pair = new Pair(ancestor0, ancestor1, affinity);
            heap.Push(pair);
          }
        }
      }
    }
    else {
      if (0) {
        static unsigned long count = 0;
      
                                                if ((count++ % 1) == 0) {
          printf("  %15.12f : %9d %9d : %15d\n", pair->affinity,
                 proxy0->children.NEntries(), proxy1->children.NEntries(),
                 heap.NEntries());
        }
      }

      // Create merged proxy
      Proxy *proxy = new Proxy(this, proxy0, proxy1);
      proxy0->cluster_representative = proxy;
      proxy1->cluster_representative = proxy;
    }

    // Delete pair
    delete pair;
  }

  // Delete merged proxies
  CleanProxies();
  // Return success
  return 1;
}

#endif // STRUCTURE_H_

