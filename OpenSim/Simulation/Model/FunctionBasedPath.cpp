#include "FunctionBasedPath.h"
#include <vector>

using namespace std;
using namespace OpenSim;
using namespace SimTK;
using SimTK::Vec3;

// Constructor
//FunctionBasedPath::Function() :
//    GeometryPath()
//{
//    setAuthors("Joris Verhagen");
//    constructProperties();
    // create interplation object
//    std::vector<OpenSim::Coordinate const*> coords;
//    for (OpenSim::Coordinate const& c :
//         this->getModel().getComponentList<OpenSim::Coordinate>()){
//        if (c.getMotionType() != OpenSim::Coordinate::MotionType::Coupled){
//            coords.push_back(&c);
//        }
//    }
//    std::vector<OpenSim::Coordinate const*> affecting_coords;
//    Nonzero_conditions cond;
//    SimTK::State st = this->getModel().initSystem();
//    for (Coordinate const* c : coords){
//        if(coord_affects_muscle(this,*c,st,cond)){
//            affecting_coords.push_back(c);
//        }
//    }
//    std::vector<int> nPoints(affecting_coords.size(),10);
//    interp = Interpolate(this,std::move(affecting_coords),st,nPoints);
//}

double FunctionBasedPath::getLength( const SimTK::State& s) const
{
//    return interp.getLength(s);
    return 0.0;
}

double FunctionBasedPath::getLengtheningSpeed( const State &s) const
{
//    return interp.getLengtheningSpeed(s);
    return 0.0;
}
