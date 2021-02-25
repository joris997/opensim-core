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
//    computePath(s);
//    return ( getCacheVariableValue<double>(s, "length") );
    return 0.0;
}

double PointBasedPath::getLengtheningSpeed( const SimTK::State& s) const
{
//    computeLengtheningSpeed(s);
//    return getCacheVariableValue<double>(s, "speed");
    return 0.0;
}
