#include "PointBasedPath.h"

using namespace std;
using namespace OpenSim;
using namespace SimTK;
using SimTK::Vec3;

// Constructor
//PointBasedPath::PointBasedPath() :
//    GeometryPath()
//{
//    setAuthors("Joris Verhagen");
//    constructProperties();
//}

double PointBasedPath::getLength( const SimTK::State& s) const
{
    computePath(s);
    return getCacheVariableValue(s, _lengthCV);
}

void PointBasedPath::setLength(const State &s, double length) const
{
    setCacheVariableValue(s, _lengthCV, length);
}


double PointBasedPath::getLengtheningSpeed( const SimTK::State& s) const
{
    computeLengtheningSpeed(s);
    return getCacheVariableValue(s, _speedCV);
}

void PointBasedPath::setLengtheningSpeed(const State &s, double speed) const
{
    setCacheVariableValue(s, _speedCV, speed);
}

//_____________________________________________________________________________
/*
 * Compute the path's moment arms for  specified coordinate.
 *
 * @param aCoord, the coordinate
 */
double PointBasedPath::computeMomentArm(const State& s, const Coordinate& aCoord) const
{
    if (!_maSolver)
        const_cast<Self*>(this)->_maSolver.reset(new MomentArmSolver(*_model));

    return _maSolver->solve(s, aCoord,  *this);
}
