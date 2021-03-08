#ifndef OPENSIM_FUNCTIONBASED_PATH_H_
#define OPENSIM_FUNCTIONBASED_PATH_H_

#include "GeometryPath.h"
#include "PointBasedPath.h"
#include "Interpolate.h"
#include <stdlib.h>

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
    Interpolate interp;
    int identity;

public:
    // Default constructor
    FunctionBasedPath();
    // Read data constructor
    FunctionBasedPath(int id);
    // Copy from PointBasedPath constructor
    FunctionBasedPath(PointBasedPath& pbp, int id);
    FunctionBasedPath(const PointBasedPath& pbp, int id);

    double getLength( const SimTK::State& s) const override;
    void setLength( const SimTK::State& s, double length) const override;
    double getLengtheningSpeed( const SimTK::State& s) const override;
    void setLengtheningSpeed( const SimTK::State& s, double speed) const override;

    double computeMomentArm(const SimTK::State &s, const Coordinate &aCoord) const override;

    void printContent();
    void readContent();
    void setIdentity( int id){identity = id;}
    int getIdentity(){return identity;}
};
}

#endif
// might need to make setLength etc. also virtual
// or make everything virtual (just for concept)
