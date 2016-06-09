//! \file Particle.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLE_H_
#define PARTICLE_H_

#define RAND1 .4 * rand() / (float) RAND_MAX - .2

#include <vector>
#include <glm/glm.hpp>

using namespace glm;

//! Common namespace for all classes related to the particle system

namespace particles {

//! Possible different particle types

enum ParticleType { FIRE = 0, SMOKE = 1 };

class ParticleSystem;

//! Individual particle of the particle system

//! Has a configurable number of sprites (display primitives) assigned to it.

class Particle {

  //! Fried declaration to allow ParticleSystem class access to protected
  //! members

  friend class ParticleSystem;

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Does nothing.

  Particle () {}

  //! Constructor for particles without sprites

  //! \param[in] position   Initial particle position
  //! \param[in] temp       Initial particle temperature
  //! \param[in] lifetime   Particle lifetime (in timesteps)

  Particle ( const vec3 &position,
             float temp,
             int lifetime ) : pos_( position ),
                              type_( FIRE ),
                              temp_( temp ),
                              lifetime_( lifetime ) {}

  //! Destructor

  //! Does nothing.

  virtual ~Particle () {}

  // ======= //
  // Getters //
  // ======= //

  //! Get distance between this and another particle

  //! \param[in] p Particle to measure distance from
  //! \return Distance between this and given particle

  float dist ( Particle& p ) {
    float dx = pos_.x - p.pos_.x;
    float dy = pos_.y - p.pos_.y;
    float dz = pos_.z - p.pos_.z;
    return sqrt( dx * dx + dy * dy + dz * dz );
  }

  //! Get current particle position

  //! \return Current particle position as vector

  const vec3 getPos() {
    return pos_;
  }

  // ======= //
  // Setters //
  // ======= //

  //! Update particle position, including assigned sprites

  //! \param[in] d Displacement vector by which particle is to be moved

  void updatePos( const vec3& d ) {
    pos_ += d;
  }


  void setSmoke( float colorCoeff) {
    type_ = SMOKE;
  }


protected:

  // ============ //
  // Data members //
  // ============ //

  //! Current particle position
  vec3 pos_;

  //! Particle type (either FIRE or SMOKE)
  ParticleType type_;

  //! Current particle temperature
  float temp_;

  //! Current particle lifetime, decreased every timestep
  int lifetime_;

};

}

#endif /* PARTICLE_H_ */
