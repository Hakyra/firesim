//! \file Emitter.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#include "Emitter.h"

namespace particles {

Emitter::Emitter ( vec3 pos,
                   vec3 size,
                   float temp,
                   int fuel,
                   int emitThreshold,
                   float fuelConsumption,
                   float lifetimeCoeff ) :
            pos_( pos ),
            size_( size ),
            temp_( temp ),
            fuel_( fuel ),
            emitThreshold_( emitThreshold ),
            fuelConsumption_( fuelConsumption ),
            lifetimeCoeff_( lifetimeCoeff ) {}

} // namespace particles
