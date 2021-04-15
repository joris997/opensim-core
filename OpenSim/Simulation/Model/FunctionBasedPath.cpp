#include "FunctionBasedPath.h"

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/PointBasedPath.h>
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>

#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <functional>
#include <sstream>

using namespace std;
using namespace OpenSim;
using namespace SimTK;
using SimTK::Vec3;

////////////////////
// Helper structs //
////////////////////
struct Discretization{
    double begin;
    double end;
    int nPoints;
    double gridsize;
};
struct Nonzero_conditions final {
    double input_val;
    double nonzero_ma;
};

// defer an action until the destruction of this wrapper
template<typename Callback>
struct Defer final {
    Callback cb;
    Defer(Callback _cb) : cb{std::move(_cb)} {
    }
    Defer(Defer const&) = delete;
    Defer(Defer&&) noexcept = default;
    Defer& operator=(Defer const&) = delete;
    Defer& operator=(Defer&&) = delete;
    ~Defer() noexcept {
        cb();
    }
};

// Define helper functions (see .cpp file)
template<typename Callback>
Defer<Callback> defer_action(Callback cb);

bool coord_affects_muscle(
    OpenSim::PointBasedPath const& pbp,
    OpenSim::Coordinate const& c,
    SimTK::State& state,
    Nonzero_conditions& out);

//////////////////////
// Helper functions //
//////////////////////
template<typename T>
void linspace(std::vector<double> &linspaced, T start_in, T end_in, int num_in){
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);

  if (num == 0){
      return;
  }
  if (num == 1){
      linspaced.push_back(start);
      return;
  }

  double delta = (end - start) / (num - 1);
  for(int i=0; i < num-1; ++i){
      linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end);
}

class Interpolate{
    private:
        // dimension of the interpolation (e.g. two states -> 2D -> 2)
        int dimension;
        // vector containing the orthogonal discretization vectors
        // e.g. disc[0] = xrange, disc[1] = yrange etc.
        std::vector<std::vector<double>> discretization;
        std::vector<double> discSizes;
        // array containing the evaluations over the discretizations
        std::vector<double> evals;
        // vector containing index closest to discretization point
        std::vector<int> n;
        // vector containing fraction of closest index to discretization
        // point to next
        std::vector<double> u;
        // array of polynomial evaluation
        std::vector<std::vector<double>> beta;
        // vector of an index in the evalsPair
        std::vector<int> loc;

        // OPENSIM INTEGRATION
        std::vector<Discretization> dS;
        std::vector<const OpenSim::Coordinate *> coords;

    public:
        // getting
        std::vector<double> getRange(int i) const {return discretization[i];}
        int getDimension() const {return dimension;}
        std::vector<Discretization> getdS() const {return dS;}
        std::vector<double> getEvals() const {return evals;}

        // template member functions
        double getEval();

        // 0th derivative
        double getInterp(const std::vector<double>& x);
        double getInterp(const SimTK::State& s);
        double getLength(const SimTK::State& s);

        // 1st derivative
        double getInterpDer(const std::vector<double>& x,
                            int coordinate, double h=0.0001);
        double getInterpDer(const SimTK::State& s,
                            int coordinate);
        double getInterpDer(const SimTK::State& s,
                            const OpenSim::Coordinate& coordinate);
        double getLengtheningSpeed(const SimTK::State& s);

        // Default constructor
        Interpolate() = default;

        // Precomputed data constructor
        explicit Interpolate(std::vector<OpenSim::Coordinate const*> coordsIn,
                             std::vector<Discretization> dSIn,
                             std::vector<double> evalsIn);

        // Constructor for non FBP/PBP related interpolation
        explicit Interpolate(std::vector<std::vector<double>> discretizationIn,
                     std::vector<std::pair<std::vector<int>,double>> evalsPair);

        // General interface constructor
        // allows one to create an interface constructor as shown below
        // with a vector of coordinates
        explicit Interpolate(OpenSim::PointBasedPath const& pbp,
                             OpenSim::Coordinate const** cBegin,
                             OpenSim::Coordinate const** cEnd,
                             SimTK::State& st,
                             std::vector<int>& nPoints);

        // Basic constructor having a vector of coordinate pointers as input
        explicit Interpolate(OpenSim::PointBasedPath const& pbp,
                             std::vector<OpenSim::Coordinate const*> coords,
                             SimTK::State& st,
                             std::vector<int>& nPoints);
};

//#define DEBUG

//////////////////////
// Helper functions //
//////////////////////
template<typename Callback>
Defer<Callback> defer_action(Callback cb) {
    return Defer<Callback>{std::move(cb)};
}

bool coord_affects_muscle(
    OpenSim::PointBasedPath const& pbp,
    OpenSim::Coordinate const& c,
    SimTK::State& state,
    Nonzero_conditions& out) {

    bool prev_locked = c.getLocked(state);
    auto reset_locked = defer_action([&] { c.setLocked(state, prev_locked); });
    double prev_val = c.getValue(state);
    auto reset_val = defer_action([&] { c.setValue(state, prev_val); });

    c.setLocked(state, false);

    static constexpr int num_steps = 3;
    double start = c.getRangeMin();
    double end = c.getRangeMax();
    double step = (end - start) / num_steps;

    for (double v = start; v <= end; v += step) {
        c.setValue(state, v);
        double ma = pbp.computeMomentArm(state, c);
        if (std::abs(ma) > 0.001) {
            out.input_val = v;
            out.nonzero_ma = ma;
            return true;
        }
    }
    return false;
}

// General vector based interpolation constructor
Interpolate::Interpolate(std::vector<OpenSim::Coordinate const*> coordsIn,
                         std::vector<Discretization> dSIn,
                         std::vector<double> evalsIn)
    : dimension(dSIn.size()),
      evals(evalsIn),
      n(dimension,0),
      u(dimension,0),
      loc(dimension,0),
      dS(dSIn),
      coords(coordsIn)
    {
        std::vector<double> dc_;
        for (int i=0; i<dimension; i++){
            beta.push_back({0,0,0,0});

            dc_.clear();
            linspace(dc_,dS[i].begin,dS[i].end,dS[i].nPoints);
            discretization.push_back(dc_);
            discSizes.push_back(dS[i].gridsize);
        }
    }

Interpolate::Interpolate(std::vector<std::vector<double>> discretizationIn,
                     std::vector<std::pair<std::vector<int>, double>> evalsPair)
    : dimension(discretizationIn.size()),
      discretization(discretizationIn),
      n(dimension,0),
      u(dimension,0),
      loc(dimension,0)
    {
        assert(discretization.size() == evalsPair[0].first.size());
        // allow it to work with the new struct method
        for (int i=0; i<dimension; i++){
            Discretization dc;
            dc.begin = discretization[i][0];
            dc.end = discretization[i][discretization[i].size()-1];
            dc.nPoints = discretization[i].size();
            dc.gridsize = (dc.end - dc.begin)/(dc.nPoints-1);
            dS.push_back(dc);
        }

        for (int i=0; i<dimension; i++){
            beta.push_back({0,0,0,0});
            discSizes.push_back(discretization[i].size());
        }
        // I'm aware this is still stupid but it is what it is for now
        int numOfLoops = 1;
        for (int i=0; i<dimension; i++){
            numOfLoops *= discretization[i].size();
        }

        std::vector<int> cnt(dimension,0);
        for (int i=0; i<numOfLoops; i++){
            for (unsigned k=0; k<evalsPair.size(); k++){
                if (evalsPair[k].first == cnt){
                    evals.push_back(evalsPair[k].second);
                }
            }
            // update cnt values
            for (int x=dimension-1; x>=0; x--){
                if (cnt[x] != dS[x].nPoints-1){
                    cnt[x] += 1;
                    break;
                }
                if (cnt[x] == dS[x].nPoints-1){
                    for (int y=x; y<dimension; y++){
                        cnt[y] = 0;
                    }
                }
            }
        }
    }

Interpolate::Interpolate(OpenSim::PointBasedPath const& gp,
                         OpenSim::Coordinate const** cBegin,
                         OpenSim::Coordinate const** cEnd,
                         SimTK::State& st,
                         std::vector<int>& nPoints)
        : dimension(nPoints.size()),
          n(dimension,0),
          u(dimension,0),
          loc(dimension,0)
    {
        // get the size of the 'vector'. Inclusive bottom, exclusive top
        std::ptrdiff_t n = (cEnd-cBegin)+1;
        assert(n > 0);
        assert(dimension == (int)n);

        // put all coordinate pointers in a vector to later unpack an incoming
        // state to a vector of coordinate values
        for (int i=0; i<dimension; i++){
            coords.push_back(cBegin[i]);
        }

        // unlock coordinates
        for (int i=0; i<dimension; i++){
            const OpenSim::Coordinate& c = *cBegin[i];
            bool c_was_locked = c.getLocked(st);
            c.setLocked(st, false);
            auto unlock_c = defer_action(
                        [&] { c.setLocked(st, c_was_locked); });
            double c_initial_value = c.getValue(st);
            auto reset_c_val = defer_action(
                        [&] { c.setValue(st, c_initial_value); });
        }

        // make discretization objects for interpolation class instance
        Discretization dc;
        for (int i=0; i<dimension; i++){
            const OpenSim::Coordinate& c = *cBegin[i];
            dc.begin = std::max(c.getRangeMin(),-(double)SimTK_PI);
            dc.end = std::min(c.getRangeMax(),(double)SimTK_PI);
            dc.nPoints = nPoints[i];
            dc.gridsize = (dc.end-dc.begin) / (dc.nPoints-1);
            dS.push_back(dc);
        }
        // slightly extend the bound for accurate interpolation on the edges
        for (int dim=0; dim<dimension; dim++){
            dS[dim].begin   -= 1*dS[dim].gridsize;
            dS[dim].end     += 2*dS[dim].gridsize;
            dS[dim].nPoints += 3;
#ifdef DEBUG
            std::cout << "dc: " << dS[dim].begin << ", "
                      << dS[dim].end << ", "
                      << dS[dim].nPoints << ", "
                      << dS[dim].gridsize << std::endl;
#endif
        }

        // compute the total number of evaluations we need to do
        int numOfLoops = 1;
        for (int i=0; i<dimension; i++){
            numOfLoops *= dS[i].nPoints;
        }

        std::vector<int> cnt(dimension);
        std::vector<double> coordValues(dimension);
        for (int i=0; i<numOfLoops; i++){
            for (int ii=0; ii<dimension; ii++){
                coordValues[ii] = dS[ii].begin + (cnt[ii]*dS[ii].gridsize);
                const OpenSim::Coordinate& c = *cBegin[ii];
                c.setValue(st,coordValues[ii]);
            }
            evals.push_back((gp.getLength(st)));

            // update cnt values
            for (int x=dimension-1; x>=0; x--){
                if (cnt[x] != dS[x].nPoints-1){
                    cnt[x] += 1;
                    break;
                } else{
                    for (int y=x; y<dimension; y++){
                        cnt[y] = 0;
                    }
                }
            }
        }
        // just make it for using the old getInterp method
        std::vector<double> dc_;
        for (int i=0; i<dimension; i++){
            beta.push_back({0,0,0,0});

            dc_.clear();
            linspace(dc_,dS[i].begin,dS[i].end,dS[i].nPoints);
            discretization.push_back(dc_);
            discSizes.push_back(dS[i].nPoints);
        }
    }

Interpolate::Interpolate(OpenSim::PointBasedPath const& pbp,
                         std::vector<OpenSim::Coordinate const*> coords,
                         SimTK::State& st,
                         std::vector<int>& nPoints)
            : Interpolate(pbp,
                          &coords[0],
                          &coords[coords.size()-1],
                          st,
                          nPoints)
            {}

/////////////
// Methods //
/////////////
// getLength based on the state
double Interpolate::getLength(const SimTK::State& s){
    return getInterp(s);
}
// getLength with mapping from State to vector x
double Interpolate::getInterp(const SimTK::State &s){
    assert(coords.size() != 0);
    std::vector<double> x;
    for (unsigned i=0; i<coords.size(); i++){
        x.push_back(coords[i]->getValue(s));
    }
    return getInterp(x);
}
double Interpolate::getInterp(const std::vector<double> &x){
    // This is the main interpolation function
    // IN:  x, a vector of points within the considered interpolation range
    // OUT: eval, the interpolated value
    assert(x.size() == dimension);

    // get the index of the closest range value to the discretization point
    for (int i=0; i<dimension; i++){
        auto it = std::find_if(std::begin(discretization[i]),
                               std::end(discretization[i]),
                               [&](double j){return j >= x[i];});
        n[i] = std::distance(discretization[i].begin(), it)-1;
    }

    // compute remaining fraction
    for (int i=0; i<dimension; i++){
        u[i] = (x[i]-discretization[i][n[i]])/
                (discretization[i][2]-discretization[i][1]);
    }

    // compute the polynomials (already evaluated)
    for (int i=0; i<dimension; i++){
        // compute binomial coefficient
        beta[i][0] = 0.5*(u[i]-1)*(u[i]-1)*(u[i]-1)*u[i]*(2*u[i]+1);
        beta[i][1] = -0.5*(u[i] - 1)*(6*u[i]*u[i]*u[i]*u[i] -
                                      9*u[i]*u[i]*u[i] + 2*u[i] + 2);
        beta[i][2] = 0.5*u[i]*(6*u[i]*u[i]*u[i]*u[i] - 15*u[i]*u[i]*u[i] +
                               9*u[i]*u[i] + u[i] + 1);
        beta[i][3] = -0.5*(u[i] - 1)*u[i]*u[i]*u[i]*(2*u[i] - 3);
    }

    // loop over all the considered points (n-dimensional) and multiply the
    // evaluation with the weight

    // create an array containing the  lengths of each discretization direction
    int discrLoopCnt[dimension];
    for (int i=0; i<dimension; i++){
        discrLoopCnt[i] = -1;
    }

    double z = 0;
    bool breakWhile = false;
    bool allTrue;
    double Beta = 1;
    while (discrLoopCnt[0] < 3){
        Beta = 1;
        for (int i=0; i<dimension; i++){
            Beta = Beta*beta[i][discrLoopCnt[i]+1];
        }

        for (int i=0; i<dimension; i++){
            loc[i] = discrLoopCnt[i] + n[i];
        }

        z += getEval()*Beta;

        // from the back to the front, check if we're already at the maximum
        // iteration on that 'nested' for loop or else increment with 1.
        // In short, everything starts with [-1,-1,-1,...] and we keep adding
        // ones until the array of the loops becomes [ 2, 2, 2, ...]
        for (int x=dimension-1; x>=0; x--){
            if (discrLoopCnt[x] != 2){
                discrLoopCnt[x] += 1;
                break;
            }
            if (discrLoopCnt[x] == 2){
                for (int y=x; y<dimension; y++){
                    discrLoopCnt[y] = -1;
                }
            }
        }

        // checking exit conditions
        if (breakWhile){
            break;
        }
        // loop through to check whether all are at max cnt or not
        allTrue = true;
        for (int i=0; i<dimension; i++){
            if (discrLoopCnt[i] != 2){
                allTrue = false;
                break;
            }
        }
        // if all true (all are 2) set breakWhile to break on the next iteration
        if (allTrue){
            breakWhile = true;
        }
    }
    return z;
}

// Derivative methods
double Interpolate::getInterpDer(const SimTK::State &s,
                                 const OpenSim::Coordinate &aCoord){
    for (unsigned i=0; i<coords.size(); i++){
        if (coords[i]->getName() == aCoord.getName()){
            return getInterpDer(s,i);
        }
    }
    return 0.0;
}

double Interpolate::getLengtheningSpeed(const SimTK::State& s){
    double lengtheningSpeed = 0;
    for (int i=0; i<dimension; i++){
        lengtheningSpeed += getInterpDer(s,i)*coords[i]->getSpeedValue(s);
    }
    return lengtheningSpeed;
}
double Interpolate::getInterpDer(const SimTK::State &s, int coordinate){
    assert(coords.size() != 0);
    std::vector<double> x;
    for (unsigned i=0; i<coords.size(); i++){
        x.push_back(coords[i]->getValue(s));
    }
    return getInterpDer(x,coordinate);
}
double Interpolate::getInterpDer(const std::vector<double> &x,
                                 int coordinate, double h){
    assert(x.size() == dimension);
    assert(coordinate <= dimension-1);
    assert(h>0 && h < (dS[coordinate].end-x[coordinate]));
    // This is the main interpolation function for derivatives
    // IN:  x, a vector of points within the considered interpolation range
    //      coordinate, the generalized coordinate of which we take derivative
    // OUT: eval, the interpolated value
    double f1 = getInterp(x);

    std::vector<double> x2 = x;
    x2[coordinate] += h;
    double f2 = getInterp(x2);

    return (f2 - f1)/h;
}

double Interpolate::getEval(){
    // get the wrapping length evaluation given a vector 'loc' which contains,
    // in ascending dimension, the index in each dimension
    int factor = 1;
    int idx = 0;

    for (int i=0; i<dimension-1; i++){
        factor=1;
        for (int ii=i+1; ii<=dimension-1; ii++){
            factor *= discSizes[ii];
        }
        idx += loc[i]*factor;
    }
    idx += loc[loc.size()-1];
    assert(idx <= evals.size());
    return evals[idx];
}

// ----- FunctionBasedPath implementation -----

struct OpenSim::FunctionBasedPath::Impl final {
    mutable Interpolate interp;

    Impl* clone() const {
        auto p = std::unique_ptr<Impl>{new Impl{}};
        p->interp = interp;
        return p.release();
    }
};

FunctionBasedPath FunctionBasedPath::fromPointBasedPath(const Model& model, const PointBasedPath& pbp)
{
    FunctionBasedPath fbp;

    // copy relevant data from source PBP
    fbp.upd_Appearance() = pbp.get_Appearance();
    fbp.setPathPointSet(pbp.getPathPointSet());
    fbp.setPathWrapSet(pbp.getWrapSet());

    Nonzero_conditions cond;

    std::unique_ptr<OpenSim::Model> modelClone{model.clone()};
    SimTK::State& initialState = modelClone->initSystem();
    modelClone->equilibrateMuscles(initialState);
    modelClone->realizeVelocity(initialState);

    std::vector<const Coordinate *> coordsThatAffectPBP;
    for (Coordinate const& c : modelClone->getComponentList<Coordinate>()){
        if (c.getMotionType() == Coordinate::MotionType::Coupled) {
            continue;
        }

        if (!coord_affects_muscle(pbp, c, initialState, cond)) {
            continue;
        }

        coordsThatAffectPBP.push_back(&c);
        std::cout << "affecting coord: " << c.getName() << std::endl;
    }

    // create vector for number of interpolation points
    std::vector<int> nPoints(coordsThatAffectPBP.size(), 5);

    // reinitialize the default-initialized interp with the points
    fbp._impl->interp = Interpolate{
            pbp,
            std::move(coordsThatAffectPBP),
            initialState,
            nPoints
    };

    return fbp;
}

FunctionBasedPath FunctionBasedPath::fromDataFile(const std::string &path)
{
    // data file for FunctionBasedPath top-level structure:
    //
    //     - plaintext
    //     - line-oriented
    //     - 3 sections (HEADERS, DISCRETIZATIONS, and EVALUATIONS) separated
    //       by blank lines
    //
    // details:
    //
    // line                                  content
    // -------------------------------------------------------------------
    // 0                                     ID (string)
    // 1                                     <empty>
    // 2                                     DIMENSIONS (int)
    // 3                                     <empty>
    // --- end HEADERS
    // 4                                     DISCRETIZATION_1
    // 5                                     DISCRETIZATION_2
    // 6-                                    DISCRETIZATION_n
    // 3+DIMENSIONS                          DISCRETIZATION_${DIMENSIONS}
    // 3+DIMENSIONS+1                        <empty>
    // --- end DISCRETIZATIONS
    // 3+DIMENSIONS+2                        EVAL_1_1
    // 3+DIMENSIONS+3                        EVAL_1_2
    // 3+DIMENSIONS+4-                       EVAL_1_m
    // 3+DIMENSIONS+${points_in_range}       EVAL_1_${points_in_range}
    // 3+DIMENSIONS+${points_in_range}+1     EVAL_2_1
    // ...
    // ? (depends on points_in_range for each discretization)
    // ?                                     EVAL_${DIMENSIONS}_${points_in_range}
    // EOF
    //
    // where DISCRETIZATION_$n is tab-delimited (\t) line that describes
    // evenly-spaced points along the `n`th DIMENSION:
    //
    // column           content
    // ------------------------------------------------
    // 0                range_start (float)
    // 1                range_end (float)
    // 2                points_in_range (int)
    // 3                grid_size (float)
    //
    // and EVAL_$n_$m is a line that contains a single real-valued evaluation
    // of the curve at:
    //
    //     range_start + (m * ((range_end - range_start)/points_in_range))
    //
    // along DIMENSION `n`

    FunctionBasedPath fbp;
    fbp.setDataPath(path);

    std::ifstream file{path};

    if (!file) {
        std::stringstream msg;
        msg << path << ": error opening FunctionBasedPath data file";
        OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
    }

    int lineno = 0;
    std::string line;

    // ID (ignored)
    while (std::getline(file, line)) {
        ++lineno;

        if (line.empty()) {
            break;
        }
    }

    // DIMENSIONS
    int dimensions = -1;
    while (std::getline(file, line)) {
        ++lineno;

        if (line.empty()) {
            if (dimensions <= -1) {
                std::stringstream msg;
                msg << path << ": L" << lineno << ": unexpected blank line (expected dimensions)";
                OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
            }
            break;
        }

        dimensions = std::stoi(line.c_str());
    }

    // [DISCRETIZATION_1..DISCRETIZATION_n]
    std::vector<Discretization> discretizations;
    discretizations.reserve(static_cast<size_t>(dimensions));
    while (std::getline(file, line)) {
        ++lineno;

        if (line.empty()) {
            if (discretizations.size() != static_cast<size_t>(dimensions)) {
                std::stringstream msg;
                msg << path << ": L" << lineno << ": invalid number of discretizations in this file: expected: " << dimensions << " got " << discretizations.size();
                OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
            }

            break;  // end of DISCRETIZATIONs
        }

        // DISCRETIZATION_i
        Discretization d;
        std::sscanf(line.c_str(), "%lf\t%lf\t%i\t%lf", &d.begin, &d.end, &d.nPoints, &d.gridsize);
        discretizations.push_back(d);
    }

    // [EVAL_1_1..EVAL_n_m]
    std::vector<double> evals;
    while (std::getline(file, line)) {
        ++lineno;

        if (line.empty()){
            break;  // end of EVALs
        }

        evals.push_back(atof(line.c_str()));
    }

    // sanity check
    {
        size_t expectedEvals = 0;
        for (const Discretization& d : discretizations) {
            expectedEvals += static_cast<size_t>(d.nPoints);
        }

        if (expectedEvals != evals.size()) {
            std::stringstream msg;
            msg << path << ": L" << lineno << ": invalid number of function evaluations in this file: expected " << expectedEvals << " got " << evals.size();
            OPENSIM_THROW(OpenSim::Exception, std::move(msg).str());
        }
    }

    // load those values into the interpolation class
    fbp._impl->interp = Interpolate{
        std::vector<const Coordinate*>(static_cast<size_t>(dimensions)),
        std::move(discretizations),
        std::move(evals)
    };

    return fbp;
}

FunctionBasedPath::FunctionBasedPath() : _impl{new Impl{}}
{
    constructProperty_data_path("");
}
FunctionBasedPath::FunctionBasedPath(const FunctionBasedPath&) = default;
FunctionBasedPath::FunctionBasedPath(FunctionBasedPath&&) = default;
FunctionBasedPath& FunctionBasedPath::operator=(const FunctionBasedPath&) = default;
FunctionBasedPath& FunctionBasedPath::operator=(FunctionBasedPath&&) = default;
FunctionBasedPath::~FunctionBasedPath() noexcept = default;

double FunctionBasedPath::getLength(const State& s) const
{
    if (isCacheVariableValid(s, _lengthCV)) {
        return getCacheVariableValue(s, _lengthCV);
    }

    // else: compute it
    double rv = _impl->interp.getLength(s);
    setCacheVariableValue(s, _lengthCV, rv);
    return rv;
}

void FunctionBasedPath::setLength(const State &s, double length) const
{
    setCacheVariableValue(s, _lengthCV, length);
}

double FunctionBasedPath::getLengtheningSpeed(const State &s) const
{
    if (isCacheVariableValid(s, _speedCV)) {
        return getCacheVariableValue(s, _speedCV);
    }

    // else: compute it
    double rv = _impl->interp.getLengtheningSpeed(s);
    setCacheVariableValue(s, _speedCV, rv);
    return rv;
}

void FunctionBasedPath::setLengtheningSpeed(const State &s, double speed) const
{
    setCacheVariableValue(s, _speedCV, speed);
}

double FunctionBasedPath::computeMomentArm(const State& s,
                                           const Coordinate& aCoord) const
{
    return _impl->interp.getInterpDer(s,aCoord);
}

void FunctionBasedPath::printContent(std::ostream& out) const
{
    // see FunctionBasedPath::fromDataFile(const std::string &path) for format
    // spec

    // HEADERS
    out << get_data_path() << '\n'
        << '\n'
        << _impl->interp.getDimension() << '\n'
        << '\n';

    // DISCRETIZATIONS
    for (const Discretization& d : _impl->interp.getdS()) {
        out << d.begin << '\t' << d.end << '\t' << d.nPoints << '\t' << d.gridsize << '\n';
    }
    out << '\n';

    // EVALS
    for (double d : _impl->interp.getEvals()) {
        out << d << '\n';
    }

    out.flush();
}
