#include "FunctionBasedPath.h"
#include "Model.h"
#include <fstream>

using namespace std;
using namespace OpenSim;
using namespace SimTK;
using SimTK::Vec3;

// Default constructor
FunctionBasedPath::FunctionBasedPath(){
    constructProperty_identity(0);
    readContent();
}

// Constructor for reading in information from file and creating
// an interpolation object based on that file
FunctionBasedPath::FunctionBasedPath(int id)
{
    constructProperty_identity(0);

    upd_identity() = id;
    readContent();
}

// Constructor used for conversion from PBP to FBP
FunctionBasedPath::FunctionBasedPath(const Model& model,
                                     const PointBasedPath& pbp,
                                     int id)
{
    constructProperty_identity(0);

    upd_identity() = id;
    upd_Appearance() = pbp.get_Appearance();
    setPathPointSet(pbp.getPathPointSet());
    setPathWrapSet(pbp.getWrapSet());

    // Constructor for copying a PointBasedPath object and creating
    // an interpolation object
    vector<const Coordinate *> coords;
    vector<const Coordinate *> affectingCoords;
    Nonzero_conditions cond;

    Model* modelClone = model.clone();
    State& stClone = modelClone->initSystem();
    modelClone->equilibrateMuscles(stClone);
    modelClone->realizeVelocity(stClone);

    // find all coordinates in the model
    for (Coordinate const& c : modelClone->getComponentList<Coordinate>()){
        if (c.getMotionType() != Coordinate::MotionType::Coupled){
            coords.push_back(&c);
        }
    }

    // find affecting coordinates for the muscle
    for (Coordinate const* c : coords){
        if (coord_affects_muscle(pbp,*c,stClone,cond)){
            affectingCoords.push_back(c);
            std::cout << "affecting coord: " << c->getName() << std::endl;
        }
    }
    // create vector for number of interpolation points
    vector<int> nPoints(affectingCoords.size(),5);
    // create interpolation object
    interp = Interpolate(pbp,move(affectingCoords),stClone,nPoints);
}


double FunctionBasedPath::getLength(const State& s) const
{
    return interp.getLength(s);
}

void FunctionBasedPath::setLength(const State &s, double length) const
{
    setCacheVariableValue(s, _lengthCV, length);
}

double FunctionBasedPath::getLengtheningSpeed(const State &s) const
{
    return interp.getLengtheningSpeed(s);
}

void FunctionBasedPath::setLengtheningSpeed(const State &s, double speed) const
{
    setCacheVariableValue(s, _speedCV, speed);
}

double FunctionBasedPath::computeMomentArm(const State& s,
                                           const Coordinate& aCoord) const
{
    return interp.getInterpDer(s,aCoord);
}


void FunctionBasedPath::printContent(std::ofstream& printFile) const{
    printFile << get_identity() << "\n";
    printFile << "\n";
    printFile << interp.getDimension() << "\n";
    printFile << "\n";

    vector<Discretization> dS = interp.getdS();
    for (unsigned i=0; i<dS.size(); i++){
        printFile << dS[i].begin << "\t" << dS[i].end << "\t" <<
                     dS[i].nPoints << "\t" << dS[i].gridsize << "\n";
    }
    printFile << "\n";

    vector<double> evals = interp.getEvals();
    for (unsigned i=0; i<evals.size(); i++){
        printFile << evals[i] << "\n";
    }

    printFile.close();
}

void FunctionBasedPath::readContent() {
    vector<const Coordinate *> coords;
    vector<Discretization> dS;
    vector<double> evals;

    ifstream readFile;
    // might change to tsv (tab/column delimeted file)
    string filename = {"FBP"+to_string(get_identity())+".xml"};
    readFile.open(filename,ios::in);

    if (readFile){
        string sLine;
        // get identity
        while(std::getline(readFile,sLine)){
            if (sLine == ""){
                break;
            }
            assert(get_identity() == atof(sLine.c_str()));
        }
        // get dimension
        while(std::getline(readFile,sLine)){
            if (sLine == ""){
                break;
            }
            int dimension = atof(sLine.c_str());
        }
        // get interpolation discretization data
        string delimiter = "\t";
        while(std::getline(readFile,sLine)){
            if (sLine == ""){
                break;
            }
            size_t last = 0;
            size_t next = 0;
            int  cnt = 0;
            Discretization disc;
            while ((next = sLine.find(delimiter,last)) != string::npos){
                if (cnt == 0){
                    disc.begin = atof(sLine.substr(last,next-last).c_str());
                }
                if (cnt == 1){
                    disc.end = atof(sLine.substr(last,next-last).c_str());
                }
                if (cnt == 2){
                    disc.nPoints = atof(sLine.substr(last,next-last).c_str());
                }
                last = next + 1;
                cnt++;
            }
            disc.gridsize = atof(sLine.substr(last).c_str());
            dS.push_back(disc);
        }
        assert(dS.size() == dimension);

        // get interpolation data
        while(std::getline(readFile,sLine)){
            if (sLine == ""){
                break;
            }
            evals.push_back(atof(sLine.c_str()));
        }
        // NEED TO FIND THE AFFECTING COORDS OF THE CORRESPONDING MUSCLE
//        for (unsigned i=0; i<coordsNames.size(); i++){
//            for (Coordinate& co : model.updComponentList<Coordinate>()){
//                if (co.getName() == coordsNames[i]){
//                    coords.push_back(co);
//                    break;
//                }
//            }
//        }
        interp = Interpolate(coords,dS,evals);
    } else {
        cout << "Could not find file" << endl;
    }
}
