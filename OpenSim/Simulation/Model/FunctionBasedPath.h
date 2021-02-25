#ifndef OPENSIM_FUNCTIONBASED_PATH_H_
#define OPENSIM_FUNCTIONBASED_PATH_H_

#include "GeometryPath.h"
#include "Interpolate.h"

#ifdef SWIG
    #ifdef OSIMSIMULATION_API
        #undef OSIMSIMULATION_API
        #define OSIMSIMULATION_API
    #endif
#endif

namespace OpenSim {

class Coordinate;
class PointForceDirection;
class ScaleSet;
class WrapResult;
class WrapObject;

class OSIMSIMULATION_API FunctionBasedPath : public GeometryPath {
//class OSIMSIMULATION_API FunctionBasedPath : public GeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPath, GeometryPath);

//    OpenSim_DECLARE_OUTPUT(length, double,
//                           getLength, SimTK::Stage::Position);
//    OpenSim_DECLARE_OUTPUT(lengthening_speed, double,
//                           getLengtheningSpeed, SimTK::Stage::Velocity);

private:
//    Interpolate interp;

public:
//    FunctionBasedPath();

    double getLength( const SimTK::State& s) const override;
    double getLengtheningSpeed( const SimTK::State& s) const override;
};
}

#endif
// might need to make setLength etc. also virtual
// or make everything virtual (just for concept)
