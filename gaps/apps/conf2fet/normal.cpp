static int
RefineNormalWithRansac(R3Vector& normal, const R3Point& origin, R3Point *points, int npoints,
  RNScalar tolerance, int niterations = 16,
  RNScalar *max_inlier_fraction = NULL, RNScalar *avg_inlier_fraction = NULL)
{
  // Check number of points
  if (npoints < 3) return 0;
  
  // Initialize inliers
  RNArray<R3Point *> inliers;
  R3Plane plane(origin, normal);
  for (int i = 0; i < npoints; i++) {
    R3Point *point = &points[i];
    RNScalar d = R3Distance(plane, *point);
    if (d > tolerance) continue; 
    inliers.Insert(point);
  }

  // Initialize best score
  R3Vector best_normal = normal;
  RNScalar score = (RNScalar) inliers.NEntries() / (RNScalar) npoints;
  RNScalar best_score = score;
  RNScalar total_score = score;
  int nscores = 1;

  // Search for a best normal
  for (int i = 0; i < niterations; i++) {
    // Guess normal by selecting three random points
    R3Point p0 = points[(int) (RNRandomScalar() * npoints) ];
    R3Point p1 = points[(int) (RNRandomScalar() * npoints) ];
    R3Point p2 = points[(int) (RNRandomScalar() * npoints) ];
    R3Vector v1 = p1 - p0;
    R3Vector v2 = p2 - p0;
    R3Vector n = v1 % v2;
    RNLength length = n.Length();
    if (RNIsZero(length)) continue;
    n /= length;
    
    // Find inliers
    inliers.Empty();
    R3Plane plane(origin, n);
    for (int i = 0; i < npoints; i++) {
      R3Point *point = &points[i];
      RNScalar d = R3Distance(plane, *point);
      if (d > tolerance) continue; 
      inliers.Insert(point);
    }

    // Compute  score
    RNScalar score = (RNScalar) inliers.NEntries() / (RNScalar) npoints;
    total_score += score;
    nscores++;

    // Check score
    if (score <= best_score) continue;
    
    // Compute normal 
    RNScalar variances[3];
    R3Point centroid = R3Centroid(inliers);
    R3Triad triad = R3PrincipleAxes(centroid, inliers, NULL, variances);
    if (RNIsZero(variances[1])) continue;
    
    // Update best stuff
    RNScalar dot = normal.Dot(triad[2]);
    best_normal = (dot >= 0) ? triad[2] : -triad[2];
    best_score = score;
  }

  // Return stuff
  normal = best_normal;
  if (max_inlier_fraction) *max_inlier_fraction = best_score;
  if (avg_inlier_fraction) *avg_inlier_fraction = total_score / nscores;
  
  // Return success
  return 1;
}


#if 0

int
RefineNormalWithICP(R3Vector& normal, const R3Point& origin, R3Point *points, int npoints, RNScalar end_tolerance, int niterations)
{
  // Initialize everything
  RNScalar start_tolerance = 4 * end_tolerance;
  RNScalar tolerance = start_tolerance;

  // Refine plane
  for (int iter = 0; iter < niterations; iter++) {
    // Create plane
    R3Plane plane(origin, normal);
    R3Vector previous_normal = normal;

    // Find inliers
    RNArray<R3Point *> inliers;
    for (int i = 0; i < npoints; i++) {
      R3Point *point = &points[i];
      RNScalar d = R3Distance(plane, *point);
      if (d > tolerance) continue;
      inliers.Insert(point);
    }

    // Check inliers
    if (inliers.NEntries() < 3) break;
    
    // Update normal 
    RNScalar variances[3];
    R3Point centroid = R3Centroid(inliers);
    R3Triad triad = R3PrincipleAxes(centroid, inliers, NULL, variances);
    if (RNIsZero(variances[1])) break;
    RNScalar dot = normal.Dot(triad[2]);
    normal = (dot >= 0) ? triad[2] : -triad[2];

    // RNScalar angle = R3InteriorAngle(normal, previous_normal);
    // if (RNIsPositive(angle)) printf("[ %d %.2f %d %.3f ]\n", iter, tolerance, inliers.NEntries(), angle);

    // Reduce tolerance
    if (iter == niterations-1) tolerance = end_tolerance;
    else tolerance = end_tolerance + 0.5 * (tolerance - end_tolerance);
  }

  // Return success
  return 1;
}
  
#endif



