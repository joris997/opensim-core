#ifndef OPENSIM_POINTBASED_PATH_H_
#define OPENSIM_POINTBASED_PATH_H_

#include "GeometryPath.h";

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

class PointBasedPath : public GeometryPath {
//class OSIMSIMULATION_API PointBasedPath : public GeometryPath {
//    OpenSim_DECLARE_CONCRETE_OBJECT(PointBasedPath, GeometryPath);

//    OpenSim_DECLARE_OUTPUT(length, double,
//                           getLength, SimTK::Stage::Position);
//    OpenSim_DECLARE_OUTPUT(lengthening_speed, double,
//                           getLengtheningSpeed, SimTK::Stage::Velocity);

public:
//    PointBasedPath();

    double getLength( const SimTK::State& s) const override;
    double getLengtheningSpeed( const SimTK::State& s) const override;
};
}

#endif
