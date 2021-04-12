#ifndef OPENSIM_FUNCTIONBASED_PATH_H_
#define OPENSIM_FUNCTIONBASED_PATH_H_

#include "GeometryPath.h"
#include "PointBasedPath.h"
#include "Interpolate.h"
#include <stdlib.h>
#include <vector>
#include <mutex>

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

    OpenSim_DECLARE_PROPERTY(identity,int,"Identity related to printed file");
//    OpenSim_DECLARE_PROPERTY(coords,std::vector<const Coordinate *>,"Coordinates related to the GeometryPath");

private:
    mutable Interpolate interp;
//    mutable std::mutex mtx;

public:
    // Default constructor
    FunctionBasedPath();
    // Read data constructor
    FunctionBasedPath(int id);
    // Copy from PointBasedPath constructor
    FunctionBasedPath(const Model& model, const PointBasedPath& pbp, int id);

    double getLength( const SimTK::State& s) const override;
    void setLength( const SimTK::State& s, double length) const override;
    double getLengtheningSpeed( const SimTK::State& s) const override;
    void setLengtheningSpeed( const SimTK::State& s, double speed) const override;

    double computeMomentArm(const SimTK::State &s, const Coordinate &aCoord) const override;

    void printContent(std::ofstream& printFile) const;
    void readContent();
    void setIdentity( int id) {upd_identity() = id;}
    int getIdentity() const {return get_identity();}
    Interpolate getInterpolate() const {return interp;}
};
}

#endif
// might need to make setLength etc. also virtual
// or make everything virtual (just for concept)
