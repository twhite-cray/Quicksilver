#ifndef MCT_HH
#define MCT_HH

#include "portability.hh"
#include "DeclareMacro.hh"
#include "Device.hh"
#include "MC_Distance_To_Facet.hh"
#include "MC_Location.hh"
#include "MC_Nearest_Facet.hh"
#include "MonteCarlo.hh"
#include "PhysicalConstants.hh"

class MC_Particle;
class Subfacet_Adjacency;

Subfacet_Adjacency &MCT_Adjacent_Facet(const MC_Location &location, MC_Particle &mc_particle, MonteCarlo* monteCarlo);

///  Calculates the distance from the specified coordinates to the
///  input segment. This is used to track to the faces of a 3D_G
///  mesh.
static inline double MCT_Nearest_Facet_3D_G_Distance_To_Segment(double plane_tolerance,
    double facet_normal_dot_direction_cosine,
    double A, double B, double C, double D,
    const double3 &facet_coords0,
    const double3 &facet_coords1,
    const double3 &facet_coords2,
    const MC_Vector &coordinate,
    const DirectionCosine *__restrict__ const direction_cosine,
    bool allow_enter)
{
  static constexpr double boundingBox_tolerance = 1e-9;
  const double numerator = -1.0*(A * coordinate.x +
      B * coordinate.y +
      C * coordinate.z +
      D);

  /* Plane equation: numerator = -P(x,y,z) = -(Ax + By + Cz + D)
if: numerator < -1e-8*length(x,y,z)   too negative!
if: numerator < 0 && numerator^2 > ( 1e-8*length(x,y,z) )^2   too negative!
reverse inequality since squaring function is decreasing for negative inputs.
If numerator is just SLIGHTLY negative, then the particle is just outside of the face */

  // Filter out too negative distances
  if (!allow_enter && numerator < 0.0 && numerator * numerator > plane_tolerance) {
    return PhysicalConstants::_hugeDouble; }

  // we have to restrict the solution to within the triangular face
  const double distance = numerator / facet_normal_dot_direction_cosine;

  // see if the intersection point of the ray and the plane is within the triangular facet
  MC_Vector intersection_pt;
  intersection_pt.x = coordinate.x + distance * direction_cosine->alpha;
  intersection_pt.y = coordinate.y + distance * direction_cosine->beta;
  intersection_pt.z = coordinate.z + distance * direction_cosine->gamma;

  // if the point is completely below the triangle, it is not in the triangle
#define IF_POINT_BELOW_CONTINUE(axis)                                    \
  if ( facet_coords0.axis > intersection_pt.axis + boundingBox_tolerance&& \
      facet_coords1.axis > intersection_pt.axis + boundingBox_tolerance && \
      facet_coords2.axis > intersection_pt.axis + boundingBox_tolerance ) { return PhysicalConstants::_hugeDouble; }

#define IF_POINT_ABOVE_CONTINUE(axis)                                    \
  if ( facet_coords0.axis < intersection_pt.axis - boundingBox_tolerance && \
      facet_coords1.axis < intersection_pt.axis - boundingBox_tolerance && \
      facet_coords2.axis < intersection_pt.axis - boundingBox_tolerance ) { return PhysicalConstants::_hugeDouble; }

  // Is the intersection point inside the triangular facet?  Project to 2D and see.

  // A^2 + B^2 + C^2 = 1, so max(|A|,|B|,|C|) >= 1/sqrt(3) = 0.577
  // (all coefficients can't be small)
  double cross0 = 0, cross1 = 0, cross2 = 0;
  if ( C < -0.5 || C > 0.5 )
  {
    IF_POINT_BELOW_CONTINUE(x);
    IF_POINT_ABOVE_CONTINUE(x);
    IF_POINT_BELOW_CONTINUE(y);
    IF_POINT_ABOVE_CONTINUE(y);

#define AB_CROSS_AC(ax,ay,bx,by,cx,cy) ( (bx-ax)*(cy-ay) - (by-ay)*(cx-ax) )

    cross1 = AB_CROSS_AC(facet_coords0.x, facet_coords0.y,
        facet_coords1.x, facet_coords1.y,
        intersection_pt.x,  intersection_pt.y);
    cross2 = AB_CROSS_AC(facet_coords1.x, facet_coords1.y,
        facet_coords2.x, facet_coords2.y,
        intersection_pt.x,  intersection_pt.y);
    cross0 = AB_CROSS_AC(facet_coords2.x, facet_coords2.y,
        facet_coords0.x, facet_coords0.y,
        intersection_pt.x,  intersection_pt.y);

  }
  else if ( B < -0.5 || B > 0.5 )
  {
    IF_POINT_BELOW_CONTINUE(x);
    IF_POINT_ABOVE_CONTINUE(x);
    IF_POINT_BELOW_CONTINUE(z);
    IF_POINT_ABOVE_CONTINUE(z);

    cross1 = AB_CROSS_AC(facet_coords0.z, facet_coords0.x,
        facet_coords1.z, facet_coords1.x,
        intersection_pt.z,  intersection_pt.x);
    cross2 = AB_CROSS_AC(facet_coords1.z, facet_coords1.x,
        facet_coords2.z, facet_coords2.x,
        intersection_pt.z,  intersection_pt.x);
    cross0 = AB_CROSS_AC(facet_coords2.z, facet_coords2.x,
        facet_coords0.z, facet_coords0.x,
        intersection_pt.z,  intersection_pt.x);

  }
  else if ( A < -0.5 || A > 0.5 )
  {
    IF_POINT_BELOW_CONTINUE(z);
    IF_POINT_ABOVE_CONTINUE(z);
    IF_POINT_BELOW_CONTINUE(y);
    IF_POINT_ABOVE_CONTINUE(y);

    cross1 = AB_CROSS_AC(facet_coords0.y, facet_coords0.z,
        facet_coords1.y, facet_coords1.z,
        intersection_pt.y,  intersection_pt.z);
    cross2 = AB_CROSS_AC(facet_coords1.y, facet_coords1.z,
        facet_coords2.y, facet_coords2.z,
        intersection_pt.y,  intersection_pt.z);
    cross0 = AB_CROSS_AC(facet_coords2.y, facet_coords2.z,
        facet_coords0.y, facet_coords0.z,
        intersection_pt.y,  intersection_pt.z);
  }

  const double cross_tol = 1e-9 * MC_FABS(cross0 + cross1 + cross2);  // cross product tolerance

  if ( (cross0 > -cross_tol && cross1 > -cross_tol && cross2 > -cross_tol) ||
      (cross0 <  cross_tol && cross1 <  cross_tol && cross2 <  cross_tol) )
  {
    return distance;
  }
  return PhysicalConstants::_hugeDouble;
}

///  Returns a coordinate that represents the "center" of the cell.
static inline MC_Vector MCT_Cell_Position_3D_G(const MC_Domain &domain, int cell_index)
{
  const int num_points = domain.mesh._cellConnectivity[cell_index].num_points;
  MC_Vector coordinate;

  for ( int point_index = 0; point_index < num_points; point_index ++ )
  {
    int point = domain.mesh._cellConnectivity[cell_index]._point[point_index];

    coordinate.x += domain.mesh._node[point].x;
    coordinate.y += domain.mesh._node[point].y;
    coordinate.z += domain.mesh._node[point].z;
  }

  const double one_over_num_points = 1.0/((double)num_points);
  coordinate.x *= one_over_num_points;
  coordinate.y *= one_over_num_points;
  coordinate.z *= one_over_num_points;

  return coordinate;
}

static inline MC_Vector MCT_Cell_Position_3D_G(const DeviceDomain &ddomain, const int cell_index)
{
  MC_Vector coordinate;

  static constexpr int num_points = DeviceCell::numQuadPoints;

  for ( int point_index = 0; point_index < num_points; point_index ++ )
  {
    const int point = ddomain.cells[cell_index].quadPoints[point_index];

    coordinate.x += ddomain.nodes[point].x;
    coordinate.y += ddomain.nodes[point].y;
    coordinate.z += ddomain.nodes[point].z;
  }

  static constexpr double one_over_num_points = 1.0/((double)num_points);
  coordinate.x *= one_over_num_points;
  coordinate.y *= one_over_num_points;
  coordinate.z *= one_over_num_points;

  return coordinate;
}

///  Move the input particle by a small amount toward the center of the cell.
static inline void MCT_Nearest_Facet_3D_G_Move_Particle(const DeviceDomain &ddomain,
    const MC_Location &location,
    MC_Vector &coordinate, // input/output: move this coordinate
    const double move_factor)      // input: multiplication factor for move
{
  MC_Vector move_to = MCT_Cell_Position_3D_G(ddomain, location.cell);

  coordinate.x += move_factor * ( move_to.x - coordinate.x );
  coordinate.y += move_factor * ( move_to.y - coordinate.y );
  coordinate.z += move_factor * ( move_to.z - coordinate.z );
}

static inline void MCT_Nearest_Facet_Find_Nearest(MC_Particle *__restrict__ const mc_particle,
    const DeviceDomain &ddomain,
    MC_Location *__restrict__ const location,
    MC_Vector &coordinate,
    int &iteration, // input/output
    double &move_factor, // input/output
    const int num_facets_per_cell,
    MC_Nearest_Facet &nearest_facet,
    int &retry /* output */ )
{

  const int max_allowed_segments = 10000000;

  retry = 0;

  if ( mc_particle )
  {
    if ( (nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble && move_factor > 0) ||
        ( mc_particle->num_segments > max_allowed_segments && nearest_facet.distance_to_facet <= 0.0 ) )
    {
      // Could not find a solution, so move the particle towards the center of the cell
      // and try again.
      MCT_Nearest_Facet_3D_G_Move_Particle(ddomain, *location, coordinate, move_factor);

      iteration++;
      move_factor *= 2.0;

      if ( move_factor > 1.0e-2 )
        move_factor = 1.0e-2;

      static constexpr int max_iterations = 10000;

      if ( iteration == max_iterations )
      {
        qs_assert(false); // If we start hitting this assertion we can
        // come up with a better mitigation plan. - dfr
        retry = 0;

      }
      else
        retry = 1;

      // Allow the distance to the current facet
      location->facet = -1;

    }
  }
}

///  Calculates the distance from the specified coordinates to each
///  of the facets of the specified cell in a three-dimensional,
///  unstructured, hexahedral (Type 3D_G) domain, storing the minimum
///  distance and associated facet number.
static inline MC_Nearest_Facet MCT_Nearest_Facet_3D_G(
    MC_Particle *__restrict__ const mc_particle,
    const DeviceDomain &ddomain,
    MC_Location &location,
    MC_Vector &coordinate,
    const DirectionCosine *__restrict__ const direction_cosine)
{
  int                    iteration = 0;
  double                 move_factor = 0.5 * PhysicalConstants::_smallDouble;

  // Initialize some data for the unstructured, hexahedral mesh.
  static constexpr int num_facets_per_cell = DeviceCell::numFacets;

  while (true) // will break out when distance is found
  {
    // Determine the distance to each facet of the cell.
    // (1e-8 * Radius)^2
    const double plane_tolerance = 1e-16*(coordinate.x*coordinate.x +
        coordinate.y*coordinate.y +
        coordinate.z*coordinate.z);

    MC_Nearest_Facet nearest_facet;

    // largest negative distance (smallest magnitude, but negative)
    MC_Nearest_Facet nearest_negative_facet;
    nearest_negative_facet.distance_to_facet = -PhysicalConstants::_hugeDouble;

    MC_Distance_To_Facet distance_to_facet;

    for (int facet_index = 0; facet_index < num_facets_per_cell; facet_index++)
    {
      distance_to_facet.distance = PhysicalConstants::_hugeDouble;

      const double4 &dplane = ddomain.cells[location.cell].facets[facet_index].plane;

      const double facet_normal_dot_direction_cosine =
        (dplane.x * direction_cosine->alpha +
         dplane.y * direction_cosine->beta +
         dplane.z * direction_cosine->gamma);

      // Consider only those facets whose outer normals have
      // a positive dot product with the direction cosine.
      // I.e. the particle is LEAVING the cell.
      if (facet_normal_dot_direction_cosine <= 0.0) { continue; }

      /* profiling with gprof showed that putting a call to MC_Facet_Coordinates_3D_G
         slowed down the code by about 10%, so we get the facet coords "by hand." */
      const int3 &dpoint = ddomain.cells[location.cell].facets[facet_index].point;

      const double3 &nodeX = ddomain.nodes[dpoint.x];
      const double3 &nodeY = ddomain.nodes[dpoint.y];
      const double3 &nodeZ = ddomain.nodes[dpoint.z];

      const double t = MCT_Nearest_Facet_3D_G_Distance_To_Segment(
          plane_tolerance,
          facet_normal_dot_direction_cosine, dplane.x, dplane.y, dplane.z, dplane.w,
          nodeX, nodeY, nodeZ,
          coordinate, direction_cosine, false);

      distance_to_facet.distance = t;

      if ( distance_to_facet.distance > 0.0 )
      {
        if ( distance_to_facet.distance <= nearest_facet.distance_to_facet )
        {
          nearest_facet.distance_to_facet = distance_to_facet.distance;
          nearest_facet.facet             = facet_index;
        }
      }
      else // zero or negative distance
      {
        if ( distance_to_facet.distance > nearest_negative_facet.distance_to_facet )
        {
          // smallest in magnitude, but negative
          nearest_negative_facet.distance_to_facet = distance_to_facet.distance;
          nearest_negative_facet.facet             = facet_index;
        }
      }
    }

    if ( nearest_facet.distance_to_facet == PhysicalConstants::_hugeDouble )
    {
      if ( nearest_negative_facet.distance_to_facet != -PhysicalConstants::_hugeDouble )
      {
        // no positive solution, so allow a negative solution, that had really small magnitude.
        nearest_facet.distance_to_facet = nearest_negative_facet.distance_to_facet;
        nearest_facet.facet             = nearest_negative_facet.facet;
      }
    }

    int retry = 0;

    MCT_Nearest_Facet_Find_Nearest(
        mc_particle, ddomain, &location, coordinate,
        iteration, move_factor, num_facets_per_cell,
        nearest_facet,
        retry);


    if (! retry) return nearest_facet;
  } // while (true)
}  // End MCT_Nearest_Facet_3D_G

///  Calculates the nearest facet of the specified cell to the
///  specified coordinates.
///
/// \return The minimum distance and facet number.
static inline MC_Nearest_Facet MCT_Nearest_Facet(MC_Particle *__restrict__ const mc_particle,
    MC_Location &location,
    MC_Vector &coordinate,
    const DirectionCosine *__restrict__ const direction_cosine,
    const double distance_threshold,
    const double current_best_distance,
    const bool new_segment,
    const Device &device )
{
  if (location.domain < 0 || location.cell < 0)
  {
    qs_assert(false);
  }
  const DeviceDomain &ddomain = device.domains[location.domain];

  MC_Nearest_Facet nearest_facet =
    MCT_Nearest_Facet_3D_G(mc_particle, ddomain, location, coordinate, direction_cosine);

  if (nearest_facet.distance_to_facet < 0) { nearest_facet.distance_to_facet = 0; }

  if (nearest_facet.distance_to_facet >= PhysicalConstants::_hugeDouble)
  {
    qs_assert(false);
    //        MC_Warning( "Infinite distance (cell not bound) for location [Reg:%d Local Dom:%d "
    //                    "Global Dom: %d Cell:%d Fac:%d], coordinate (%g %g %g) and direction (%g %g %g).\n",
    //                    location.region, location.domain,
    //                    mcco->region->Global_Domain_Number(location.region, location.domain),
    //                    location.cell, location.facet,
    //                    coordinate.x, coordinate.y, coordinate.z,
    //                    direction_cosine->alpha, direction_cosine->beta, direction_cosine->gamma);
    //        if ( mc_particle )
    //        {
    //           MC_Warning( "mc_particle.identifier %" PRIu64 "\n", mc_particle->identifier );
    //        }
  }

  return nearest_facet;
}  // End MCT_Nearest_Facet

///  Fills in the facet_points array with the domain local point
///  numbers specified by the cell number and cell-local facet number
///  for a 3DG mesh.
static inline void MCT_Facet_Points_3D_G(const MC_Domain    &domain,               // input
    const int                 cell,                 // input
    const int                 facet,                // input
    const int                 num_points_per_facet, // input
    int                *__restrict__ const facet_points          /* output */)
{
  // Determine the domain local points of the facet in the cell for the 2DG or 3DG mesh.
  for ( int point_index = 0; point_index < num_points_per_facet; point_index++ )
    facet_points[point_index] = domain.mesh._cellConnectivity[cell]._facet[facet].point[point_index];
}

///  \return 6 times the volume of the tet.
///
///  subtract v3 from v0, v1 and v2.  Then take the triple product of v0, v1 and v2.
static inline double MCT_Cell_Volume_3D_G_vector_tetDet(const MC_Vector &v0_,
    const MC_Vector &v1_,
    const MC_Vector &v2_,
    const MC_Vector &v3)
{
  MC_Vector v0(v0_), v1(v1_), v2(v2_);

  v0.x -= v3.x; v0.y -= v3.y; v0.z -= v3.z;
  v1.x -= v3.x; v1.y -= v3.y; v1.z -= v3.z;
  v2.x -= v3.x; v2.y -= v3.y; v2.z -= v3.z;

  return
    v0.z*(v1.x*v2.y - v1.y*v2.x) +
    v0.y*(v1.z*v2.x - v1.x*v2.z) +
    v0.x*(v1.y*v2.z - v1.z*v2.y);
}

///  Generates a random coordinate inside a polyhedral cell.
static inline void MCT_Generate_Coordinate_3D_G(uint64_t *random_number_seed,
    int domain_num,
    int cell,
    MC_Vector &coordinate,  
    MonteCarlo* monteCarlo )
{
  const MC_Domain &domain = monteCarlo->domain[domain_num];

  // Determine the cell-center nodal point coordinates.
  MC_Vector center = MCT_Cell_Position_3D_G(domain, cell);

  int num_facets = domain.mesh._cellConnectivity[cell].num_facets;
  if (num_facets == 0)
  {
    coordinate.x = coordinate.y = coordinate.z = 0;
    return;
  }

  double random_number = rngSample(random_number_seed);
  double which_volume = random_number * 6.0 * domain.cell_state[cell]._volume;

  // Find the tet to sample from.
  double current_volume = 0.0;
  int facet_index = -1;
  const MC_Vector *point0 = NULL;
  const MC_Vector *point1 = NULL;
  const MC_Vector *point2 = NULL;
  while (current_volume < which_volume)
  {
    facet_index++;

    if (facet_index == num_facets) { break; }

    int facet_points[3];
    MCT_Facet_Points_3D_G(domain, cell, facet_index, 3, facet_points);
    point0 = &domain.mesh._node[facet_points[0]];
    point1 = &domain.mesh._node[facet_points[1]];
    point2 = &domain.mesh._node[facet_points[2]];

    double subvolume = MCT_Cell_Volume_3D_G_vector_tetDet(*point0, *point1, *point2, center);
    current_volume += subvolume;

  }

  // Sample from the tet.
  double r1 = rngSample(random_number_seed);
  double r2 = rngSample(random_number_seed);
  double r3 = rngSample(random_number_seed);

  // Cut and fold cube into prism.
  if (r1 + r2 > 1.0)
  {
    r1 = 1.0 - r1;
    r2 = 1.0 - r2;
  }
  // Cut and fold prism into tetrahedron.
  if (r2 + r3 > 1.0)
  {
    double tmp = r3;
    r3 = 1.0 - r1 - r2;
    r2 = 1.0 - tmp;
  }
  else if (r1 + r2 + r3 > 1.0)
  {
    double tmp = r3;
    r3 = r1 + r2 + r3 - 1.0;
    r1 = 1.0 - r2 - tmp;
  }

  // numbers 1-4 are the barycentric coordinates of the random point.
  double r4 = 1.0 - r1 - r2 - r3;

  // error check
  if ((point0 == NULL) || (point1 == NULL) || (point2 == NULL))
  {
    MC_Fatal_Jump( "Programmer Error: points must not be NULL: point0=%p point1=%p point2=%p",
        point0, point1, point2);
    return;
  }

  coordinate.x = ( r4 * center.x + r1 * point0->x + r2 * point1->x + r3 * point2->x );
  coordinate.y = ( r4 * center.y + r1 * point0->y + r2 * point1->y + r3 * point2->y );
  coordinate.z = ( r4 * center.z + r1 * point0->z + r2 * point1->z + r3 * point2->z );
}

///  Reflects the particle off of a reflection boundary.
static inline void MCT_Reflect_Particle(const Device &device, MC_Particle &particle)
{
  DirectionCosine *__restrict__ const direction_cosine = particle.Get_Direction_Cosine();
  MC_Location              location = particle.Get_Location();

  const double4 &facet_normal = device.domains[location.domain].cells[location.cell].facets[location.facet].plane;

  const double dot = 2.0*( direction_cosine->alpha * facet_normal.x +
      direction_cosine->beta  * facet_normal.y +
      direction_cosine->gamma * facet_normal.z );

  if ( dot > 0 ) // do not reflect a particle that is ALREADY pointing inward
  {
    // reflect the particle
    direction_cosine->alpha -= dot * facet_normal.x;
    direction_cosine->beta  -= dot * facet_normal.y;
    direction_cosine->gamma -= dot * facet_normal.z;
  }

  // Calculate the reflected, velocity components.
  const double particle_speed = particle.velocity.Length();
  particle.velocity.x = particle_speed * particle.direction_cosine.alpha;
  particle.velocity.y = particle_speed * particle.direction_cosine.beta;
  particle.velocity.z = particle_speed * particle.direction_cosine.gamma;
}

#endif
