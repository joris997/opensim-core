#include "FunctionBasedPath.h"
#include "Model.h"
#include <vector>
#include <fstream>

using namespace std;
using namespace OpenSim;
using namespace SimTK;
using SimTK::Vec3;

FunctionBasedPath::FunctionBasedPath(){

}

FunctionBasedPath::FunctionBasedPath(const int id){
    // Constructor for reading in information from file and creating
    // an interpolation object based on that file
    identity = id;
    readContent();
}

FunctionBasedPath::FunctionBasedPath(PointBasedPath& pbp, int id){
    identity = id;
    // Constructor for copying a PointBasedPath object and creating
    // an interpolation object
    vector<Coordinate const*> coords;
    vector<Coordinate const*> affectingCoords;
    Nonzero_conditions cond;

    // find all coordinates in the model
    const Model modelConst = pbp.getModel();
    Model model = modelConst;
    for (Coordinate const& c : model.getComponentList<Coordinate>()){
        if (c.getMotionType() != Coordinate::MotionType::Coupled){
            coords.push_back(&c);
        }
    }
    model.buildSystem();
    State& st = model.initializeState();

    // find affecting coordinates for the muscle
    for (Coordinate const* c : coords){
        if (coord_affects_muscle(pbp,*c,st,cond)){
            affectingCoords.push_back(c);
        }
    }
    // create vector for number of interpolation points
    vector<int> nPoints(affectingCoords.size(),10);
    // create interpolation object
    interp = Interpolate(pbp,move(affectingCoords),st,nPoints);
}

FunctionBasedPath::FunctionBasedPath(const PointBasedPath& pbp, int id){
    identity = id;
    // Constructor for copying a PointBasedPath object and creating
    // an interpolation object
    vector<Coordinate const*> coords;
    vector<Coordinate const*> affectingCoords;
    Nonzero_conditions cond;

    // find all coordinates in the model
    const Model modelConst = pbp.getModel();
    Model model = modelConst;
    for (Coordinate const& c : model.getComponentList<Coordinate>()){
        if (c.getMotionType() != Coordinate::MotionType::Coupled){
            coords.push_back(&c);
        }
    }
    model.buildSystem();
    State& st = model.initializeState();

    // find affecting coordinates for the muscle
    for (Coordinate const* c : coords){
        if (coord_affects_muscle(pbp,*c,st,cond)){
            affectingCoords.push_back(c);
        }
    }
    // create vector for number of interpolation points
    vector<int> nPoints(affectingCoords.size(),10);
    // create interpolation object
    interp = Interpolate(pbp,move(affectingCoords),st,nPoints);
}

double FunctionBasedPath::getLength(const State& s) const
{
    Interpolate interpCopy = interp;
    return interpCopy.getLength(s);
}

void FunctionBasedPath::setLength(const State &s, double length) const
{
    setCacheVariableValue(s, _lengthCV, length);
}

double FunctionBasedPath::getLengtheningSpeed(const State &s) const
{
    Interpolate interpCopy = interp;
    return interpCopy.getLengtheningSpeed(s);
}

void FunctionBasedPath::setLengtheningSpeed(const State &s, double speed) const
{
    setCacheVariableValue(s, _speedCV, speed);
}

double FunctionBasedPath::computeMomentArm(const State& s, const Coordinate& aCoord) const
{
    if (!_maSolver)
        const_cast<Self*>(this)->_maSolver.reset(new MomentArmSolver(*_model));

    return _maSolver->solve(s, aCoord,  *this);
}


void FunctionBasedPath::printContent(){
    ofstream printFile;
    string filename = {"FBP"+to_string(identity)+".xml"};
    printFile.open(filename);
    printFile << identity << "\n";
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

void FunctionBasedPath::readContent(){
    int dimension;
    vector<Coordinate const*> coords;
    vector<Discretization> dS;
    vector<double> evals;

    ifstream readFile;
    string filename = {"FBP"+to_string(identity)+".xml"};
    readFile.open(filename,ios::in);

    if (readFile){
        string sLine;
        // get identity
        while(std::getline(readFile,sLine)){
            if (sLine == ""){
                break;
            }
            assert(identity == atof(sLine.c_str()));
        }
        // get dimension
        while(std::getline(readFile,sLine)){
            if (sLine == ""){
                break;
            }
            dimension = atof(sLine.c_str());
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
        interp = Interpolate(coords,dS,evals);
    } else {
        cout << "Could not find file" << endl;
    }
}
