#ifndef MC_PARTICLE_INCLUDE
#define MC_PARTICLE_INCLUDE

#include <cinttypes>

#include "DirectionCosine.hh"
#include "MC_Tally_Event.hh"

#include "MC_Vector.hh"
#include "MC_Facet_Adjacency.hh"
#include "MC_Location.hh"

#include "DeclareMacro.hh"

class MC_Base_Particle;


class MC_Particle
{
 public:

    // the current position of the particle
    MC_Vector coordinate;

    // the velocity of the particle
    MC_Vector velocity;

    // the direction of the particle
    DirectionCosine direction_cosine;

    // the kinetic energy of the particle
    double kinetic_energy;

    // the weight of the particle
    double weight;

    // the time remaining for this particle to hit census
    double time_to_census;

    // cacheing the current total cross section
    double totalCrossSection;

    // the age of this particle
    double age;

    // the number of mean free paths to a collision
    double num_mean_free_paths;

    // distance to a collision
    double mean_free_path;

    // the distance this particle travels in a segment.
    double segment_path_length;

    // the random number seed for the random number generator for this particle
    uint64_t random_number_seed;

    // unique identifier used to identify and track individual particles in the simulation
    uint64_t identifier;

   // the last event this particle underwent
    MC_Tally_Event::Enum last_event;

    int num_collisions;

    double num_segments;

    // the breed of the particle how it was produced
    int breed;

    // current energy group of the particle
    int energy_group;

    // its current domain in the spatial decomposition
    int domain;

    // the current cell in its current domain
    int cell;

    int facet;

    // When crossing a facet, keep the surface normal dot product
    double normal_dot;

public:
   __host__ __device__ MC_Particle();

   MC_Particle( const MC_Base_Particle &from_particle );

   bool operator==(const MC_Particle &that) const
   {
     assert(coordinate == that.coordinate);
     assert(velocity == that.velocity);
     assert(direction_cosine == that.direction_cosine);
     assert(kinetic_energy == that.kinetic_energy);
     assert(weight == that.weight);
     assert(time_to_census == that.time_to_census);
     assert(totalCrossSection == that.totalCrossSection);
     assert(age == that.age);
     assert(num_mean_free_paths == that.num_mean_free_paths);
     assert(mean_free_path == that.mean_free_path);
     assert(segment_path_length == that.segment_path_length);
     assert(random_number_seed == that.random_number_seed);
     assert(identifier == that.identifier);
     assert(last_event == that.last_event);
     assert(num_collisions == that.num_collisions);
     assert(num_segments == that.num_segments);
     assert(breed == that.breed);
     assert(energy_group == that.energy_group);
     assert(domain == that.domain);
     assert(cell == that.cell);
     assert(facet == that.facet);
     assert(normal_dot == that.normal_dot);
      return (
         (coordinate == that.coordinate) &&
         (velocity == that.velocity) &&
         (direction_cosine == that.direction_cosine) &&
         (kinetic_energy == that.kinetic_energy) &&
         (weight == that.weight) &&
         (time_to_census == that.time_to_census) &&
         (totalCrossSection == that.totalCrossSection) &&
         (age == that.age) &&
         (num_mean_free_paths == that.num_mean_free_paths) &&
         (mean_free_path == that.mean_free_path) &&
         (segment_path_length == that.segment_path_length) &&
         (random_number_seed == that.random_number_seed) &&
         (identifier == that.identifier) &&
         (last_event == that.last_event) &&
         (num_collisions == that.num_collisions) &&
         (num_segments == that.num_segments) &&
         (breed == that.breed) &&
         (energy_group == that.energy_group) &&
         (domain == that.domain) &&
         (cell == that.cell) &&
         (facet == that.facet) &&
         (normal_dot == that.normal_dot));
   }

   void Copy_From_Base( const MC_Base_Particle &from_particle);

   __host__ __device__ MC_Location Get_Location() const;

   // format a string with the contents of the particle
   void Copy_Particle_To_String(std::string &output_string) const;

   // move a particle a distance in the direction_cosine direction
   __host__ __device__ void Move_Particle(const DirectionCosine & direction_cosine, const double distance);

   void PrintParticle();
};

//----------------------------------------------------------------------------------------------------------------------
//  Return a MC_Location given domain, cell, facet.
//----------------------------------------------------------------------------------------------------------------------
__host__ __device__ inline MC_Location MC_Particle::Get_Location() const
{
    return MC_Location(domain, cell, facet);
}

//----------------------------------------------------------------------------------------------------------------------
//  Move the particle a straight-line distance along a specified cosine.
//----------------------------------------------------------------------------------------------------------------------
__host__ __device__ inline void MC_Particle::Move_Particle( const DirectionCosine &my_direction_cosine,
                                      const double distance)
{
    coordinate.x += (my_direction_cosine.alpha * distance);
    coordinate.y += (my_direction_cosine.beta  * distance);
    coordinate.z += (my_direction_cosine.gamma * distance);
}

//----------------------------------------------------------------------------------------------------------------------
//  Print all of the particles components
//----------------------------------------------------------------------------------------------------------------------
inline void MC_Particle::PrintParticle()
{
    printf( "coordiante:          %g\t%g\t%g\n", coordinate.x, coordinate.y, coordinate.z );
    printf( "velocity:            %g\t%g\t%g\n", velocity.x, velocity.y, velocity.z );
    printf( "direction_cosine:    %g\t%g\t%g\n", direction_cosine.alpha, direction_cosine.beta, direction_cosine.gamma );
    printf( "kinetic_energy:      %g\n", kinetic_energy );
    printf( "Weight:              %g\n", weight);
    printf( "time_to_census:      %g\n", time_to_census);
    printf( "totalCrossSection:   %g\n", totalCrossSection);
    printf( "age:                 %g\n", age);
    printf( "num_mean_free_paths: %g\n", num_mean_free_paths);
    printf( "mean_free_path:      %g\n", mean_free_path);
    printf( "segment_path_length: %g\n", segment_path_length);
    printf( "random_number_seed:  %" PRIu64 "\n", random_number_seed);
    printf( "identifier:          %" PRIu64 "\n", identifier);
    printf( "last_event:          %d\n", last_event);
    printf( "num_collision:       %d\n", num_collisions);
    printf( "num_segments:        %g\n", num_segments);
    printf( "breed:               %d\n", breed);
    printf( "energy_group:        %d\n", energy_group);
    printf( "domain:              %d\n", domain);
    printf( "cell:                %d\n", cell);
    printf( "facet:               %d\n", facet);
    printf( "normal_dot:          %g\n", normal_dot);
    printf("\n");
}


int                         MC_Copy_Particle_Get_Num_Fields();


#endif //  MC_PARTICLE_INCLUDE
