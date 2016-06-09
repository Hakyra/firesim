//! \file ParticleSystem.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLESYSTEM_H_
#define PARTICLESYSTEM_H_

#ifndef RealType
#define RealType float
#endif

#include <SDL2/SDL.h>

#include "../lbm/LBM.h"
#include "Emitter.h"
#include <glm/glm.hpp>

//! Common namespace for all classes related to the particle system

namespace particles {

struct Sphere {

  Sphere( float x,
          float y,
          float z,
          float radius,
          float u_x = 0.0,
          float u_y = 0.0,
          float u_z = 0.0) :
      pos( x, y, z ), r( radius ), u( u_x, u_y, u_z ) {}

  void move() {
    pos += u;
  }

  bool isPointInside( const vec3& p ) {
    return   (p.x - pos.x) * (p.x - pos.x)
           + (p.y - pos.y) * (p.y - pos.y)
           + (p.z - pos.z) * (p.z - pos.z) < r * r;
  }

    vec3 pos;
  float r;
    vec3 u;
};

//! Particle system that handles creation, movement and visualization of
//! particles.

class ParticleSystem {

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Does nothing.

  ParticleSystem() {}

  //! Constructor that initializes the particle system according to given
  //! configuration file

  //! \param[in] configFileName Path to configuration file to parse

  ParticleSystem ( std::string configFileName );

  //! Destructor

  //! Does nothing.

  virtual ~ParticleSystem () {}

  //! Set up the particle system according to configuration options

  //! \param[in] base Root block of the parsed configuration file

  void setup( ConfBlock& base );

  //! Main simulation loop

  void run();


    void InitSDL();

    void Draw();

    void CloseSDL();


protected:


    SDL_Window* fenetre;
    SDL_GLContext contexteOpenGL;


  // ========================= //
  // Internal helper functions //
  // ========================= //

  //! Go over all particles and update position, color and size

  inline void updateParticles();

  //! Go over all emitters and emit new particles

  inline void emitParticles();

  //! Generate a black body color table

  //! \param[in] maxTemp Maximum emitted particle temperature

  void generateBlackBodyColorTable( float maxTemp );

  float getTime( timeval &start, timeval &end ) {
    return (float) ( end.tv_sec - start.tv_sec )
            + (float) ( end.tv_usec - start.tv_usec ) / 1000000.;
  }

  // ============ //
  // Data members //
  // ============ //

  //! Lattice Boltzmann fluid solver
  lbm::LBM<float> solver_;

  //! Vector containing all emitters of the particle system
  std::vector< Emitter > emitters_;

  //! Domain size in x-direction
  int sizeX_;

  //! Domain size in y-direction
  int sizeY_;

  //! Domain size in z-direction
  int sizeZ_;

  //! Total simulation timesteps
  int maxSteps_;

  //! Total number of generated particles
  int numParticles_;

  //! Thermal conservation coefficient
  float alpha_;

  //! Thermal transferability coefficient
  float beta_;

  //! Temperature threshold for a fire particle to turn into a smoke particle
  float smokeTemp_;

  //! Ambient temperature
  float ambTemp_;

  //! Thermal expansion coefficient
  float k_;

  //! Vector of inversed gravity
  vec3 gravity_;

  //! Basic particle size
  float sizeBase_;

  //! Variable particle size, multiplied with a factor depending on lifetime
  float sizeVar_;

  //! Number of sprites assigned to each particle
  int numSprites_;

  //! Table of a precomputed Gaussian distribution
  std::vector<float> gaussTable_;

  std::vector< Sphere > spheres_;

  std::string updFileName_;

    vec3 minPart;
    vec3 maxPart;

  bool dynamicLights_;
};

} // namespace particles

#endif /* PARTICLESYSTEM_H_ */
