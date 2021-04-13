#include "Interpolate.h"

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

//////////////////
// Constructors //
//////////////////
// Default constructor
Interpolate::Interpolate(){};

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
    : discretization(discretizationIn),
      dimension(discretizationIn.size()),
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
                         std::vector<int>& discretizationNPoints)
        : dimension(discretizationNPoints.size()),
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
            dc.nPoints = discretizationNPoints[i];
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
                         std::vector<int>& discretizationNPoints)
            : Interpolate(pbp,
                          &coords[0],
                          &coords[coords.size()-1],
                          st,
                          discretizationNPoints)
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
