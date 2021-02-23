#include "Interpolate.h"

template<typename T>
void printVector(std::vector<T> x){
    for (int i=0; i<x.size(); i++){
        std::cout << "x[" << i << "]: " << x[i] << std::endl;
    }
}


int binomialCoefficient(int n, int k) {
    if (k == 0 || k == n){
        return 1;
    }
    return binomialCoefficient(n-1, k-1) + binomialCoefficient(n-1, k);
}

////////////////////
// FUNCTIONS ADAM //
////////////////////
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
// factory function to Defer
template<typename Callback>
Defer<Callback> defer_action(Callback cb) {
    return Defer<Callback>{std::move(cb)};
}
static bool coord_affects_muscle(
    OpenSim::Muscle const& m,
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
        double ma = m.getGeometryPath().computeMomentArm(state, c);
        if (std::abs(ma) > 0.001) {
            out.input_val = v;
            out.nonzero_ma = ma;
            return true;
        }
    }
    return false;
}

//////////////////
// CONSTRUCTORS //
//////////////////
Interpolate::Interpolate(std::vector<std::vector<double> > discretizationIn,
                         std::vector<std::pair<std::vector<int>, double> > evalsPair)
    : discretization(discretizationIn),
      dimension(discretizationIn.size()),
      evalsPair(evalsPair),
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
            for (int k=0; k<evalsPair.size(); k++){
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

Interpolate::Interpolate(OpenSim::Muscle const& m,
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

        OpenSim::GeometryPath const& musc_path = m.getGeometryPath();

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
            auto unlock_c = defer_action([&] { c.setLocked(st, c_was_locked); });
            double c_initial_value = c.getValue(st);
            auto reset_c_val = defer_action([&] { c.setValue(st, c_initial_value); });
        }

        // make discretization objects for interpolation class instance
        Discretization dc;
        for (int i=0; i<dimension; i++){
            const OpenSim::Coordinate& c = *cBegin[i];
            dc.begin = c.getRangeMin();
            dc.end = c.getRangeMax();
            dc.nPoints = discretizationNPoints[i];
            dc.gridsize = (dc.end-dc.begin) / dc.nPoints;
            dS.push_back(dc);
        }
        // slightly extend the bound for accurate interpolation on the edges
        for (int dim=0; dim<dimension; dim++){
            dS[dim].begin   -= 1*dS[dim].gridsize;
            dS[dim].end     += 2*dS[dim].gridsize;
            dS[dim].nPoints += 3;
        }


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
            evals.push_back((musc_path.getLength(st)));
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

Interpolate::Interpolate(OpenSim::Muscle const& m,
                         std::vector<OpenSim::Coordinate const*> coords,
                         SimTK::State& st,
                         std::vector<int>& discretizationNPoints)
            : Interpolate(m,
                          &coords[0],
                          &coords[coords.size()-1],
                          st,
                          discretizationNPoints)
            {}

/////////////////
// INTERPOLATE //
/////////////////

double Interpolate::getInterp(const SimTK::State &s){
    assert(coords.size() != 0);
    std::vector<double> x;
    for (int i=0; i<coords.size(); i++){
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
        beta[i][0] = (0.5*pow(u[i] - 1,3)*u[i]*(2*u[i] + 1));
        beta[i][1] = (-0.5*(u[i] - 1)*(6*pow(u[i],4) - 9*pow(u[i],3) + 2*u[i] + 2));
        beta[i][2] = (0.5*u[i]*(6*pow(u[i],4) - 15*pow(u[i],3) + 9*pow(u[i],2) + u[i] + 1));
        beta[i][3] = (-0.5*(u[i] - 1)*pow(u[i],3)*(2*u[i] - 3));
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
//    int cnt=0;
//    while (cnt < pow(dimension,4)){
    while (discrLoopCnt[0] < 3){
        Beta = 1;
        for (int i=0; i<dimension; i++){
            Beta = Beta*beta[i][discrLoopCnt[i]+g];
        }

        for (int i=0; i<dimension; i++){
            loc[i] = discrLoopCnt[i] + n[i];
        }

        z += getEvalFast()*Beta;


        // from the back to the front, check if we're already at the maximum iteration
        // on that 'nested' for loop or else increment with 1. In short, everything
        // starts with [-1,-1,-1,...] and we keep adding ones until the array of the
        // loops becomes [ 2, 2, 2, ...]
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

        // CHECKING EXIT CONDITIONS
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
        // if all are true (all are 2) set breakWhile to break on the next iteration
        if (allTrue){
            breakWhile = true;
        }
    }
    return z;
}


double Interpolate::getInterpStruct(const SimTK::State &s){
    assert(coords.size() != 0);
    std::vector<double> x;
    for (int i=0; i<coords.size(); i++){
        x.push_back(coords[i]->getValue(s));
    }
    return getInterpStruct(x);
}
double Interpolate::getInterpStruct(const std::vector<double> &x){
    assert(x.size() == dimension);

    // get the index of the closest range value to the discretization point
    for (int i=0; i<dimension; i++){
        n[i] = floor((x[i]-dS[i].begin)/dS[i].gridsize);
    }

    // compute remaining fraction
    for (int i=0; i<dimension; i++){
        u[i] = (x[i]-(dS[i].begin + n[i]*dS[i].gridsize))/(dS[i].gridsize);
    }

    // compute the polynomials (already evaluated)
    for (int i=0; i<dimension; i++){
        // compute binomial coefficient
        beta[i][0] = (0.5*pow(u[i] - 1,3)*u[i]*(2*u[i] + 1));
        beta[i][1] = (-0.5*(u[i] - 1)*(6*pow(u[i],4) - 9*pow(u[i],3) + 2*u[i] + 2));
        beta[i][2] = (0.5*u[i]*(6*pow(u[i],4) - 15*pow(u[i],3) + 9*pow(u[i],2) + u[i] + 1));
        beta[i][3] = (-0.5*(u[i] - 1)*pow(u[i],3)*(2*u[i] - 3));
    }

    int discrLoopCnt[dimension];
    for (int i=0; i<dimension; i++){
        discrLoopCnt[i] = -1;
    }

    double z = 0;
    double Beta = 1;
    for (int cnt=0; cnt<pow(4,dimension); cnt++){
        Beta = 1;
        for (int i=0; i<dimension; i++){
            Beta *= beta[i][discrLoopCnt[i]+g];
        }

        for (int i=0; i<dimension; i++){
            loc[i] = discrLoopCnt[i] + n[i];
        }

        z += getEvalFast()*Beta;

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
    }
    return z;
}



double Interpolate::getLengtheningSpeed(const SimTK::State& s){
    double lengtheningSpeed = 0;
    for (int i=0; i<dimension; i++){
        lengtheningSpeed += getInterpDer(s,i);
    }
    return lengtheningSpeed;
}
double Interpolate::getInterpDer(const SimTK::State &s, int coordinate, double h){
    assert(coords.size() != 0);
    std::vector<double> x;
    for (int i=0; i<coords.size(); i++){
        x.push_back(coords[i]->getValue(s));
    }
    return getInterpDer(x,coordinate,h);
}
double Interpolate::getInterpDer(const std::vector<double> &x, int coordinate, double h){
    assert(x.size() == dimension);
    assert(coordinate <= dimension-1);
    assert(h>0 && h < (dS[coordinate].end-x[coordinate]));
    // This is the main interpolation function
    // IN:  x, a vector of points within the considered interpolation range
    //      coordinate, the generalized coordinate of which we take the derivative
    // OUT: eval, the interpolated value
    double f1 = getInterp(x);

    std::vector<double> x2 = x;
    x2[coordinate] += h;
    double f2 = getInterp(x2);

    return (f2 - f1)/h;
}



double Interpolate::getEvalFast(){
    // get the wrapping length evaluation given a vector 'loc' which contains, in
    // ascending dimension, the index in each dimension
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


double Interpolate::getEval(){
    // get the wrapping length evaluation given a vector 'loc' which contains, in
    // ascending dimension, the index in each dimension
    for (int i=0; i<evalsPair.size(); i++){
        if (evalsPair[i].first == loc){
            return evalsPair[i].second;
        }
    }
    return 0.0;
}


void Interpolate::computeBasisFunctions(std::vector<double> u,
                                   int order){
    // create the basis spline functions
    for (int i=0; i<dimension; i++){
        double Bk[4];
        for (int k=0; k<4; k++){
            Bk[k] = binomialCoefficient(order,k)*pow(u[i],k)*pow(1-u[i],order-k);
        }
        std::vector<double> betaArray;
        betaArray.push_back(Bk[0] + Bk[1]);
        betaArray.push_back(Bk[1]/3);
        betaArray.push_back(Bk[3] + Bk[2]);
        betaArray.push_back(-Bk[2]/3);
        beta.push_back(betaArray);
    }
}



//double z, Beta;
//int discrLoopCnt = -1;
//for (int i=0; i<4; i++){
//    Beta = beta[coordinate][i];
//    for (int j=0; j<dimension; j++){
//        loc[i] = n[j];
//    }
//    loc[coordinate] += discrLoopCnt;

//    z += getEval()*Beta/
//            (discretization[coordinate][1]-discretization[coordinate][0]);

//    discrLoopCnt += 1;
//}

//double Interpolate::interpCubicHermiteSpline(vector<double> x, int derivativeOrder){
//    // This is the main interpolation function
//    // IN:  x, a vector of points within the considered interpolation range
//    // OUT: eval, the interpolated value
//    vector<int> n, pk, pkp1, pkm1, pkp2;
//    vector<double> u, d;
//    vector<vector<double>> beta, m;
//    double Beta0=1, Beta1=1, Beta2=1, Beta3=1, z=0;
//    int order = 3;

//    // check size of input argument
//    if (x.size() != dimension){
//        std::cout << "Interpolation state not equal to object's dimension" << std::endl;
//        exit(EXIT_FAILURE);
//    }

//    // get the index of the closest range value to the discretization point
//    for (int i=0; i<dimension; i++){
//        for (int ii=1; ii<discretization[i].size(); ii++){
//            if (x[i] - discretization[i][ii] < 0){
//                n.push_back(ii-1);
//                break;
//            }
//        }
//    }
////    for (int i=0; i<n.size(); i++){
////        std::cout << "n[" << i << "]: " << n[i] << std::endl;
////    }

//    // compute remaining fraction
//    for (int i=0; i<dimension; i++){
//        d.push_back((double)(discretization[i][2]-discretization[i][1]));
//    }
//    for (int i=0; i<dimension; i++){
//        u.push_back((double)((x[i]-discretization[i][n[i]])/(d[i])));
//    }

//    // compute the polynomials (already evaluated)
//    if (derivativeOrder == 0){
//        computeBasisFunctions(beta,u,order);
//    } else if (derivativeOrder == 1){
//        computeBasisFunctionsDerivatives(beta,u,order);
//    } else {
//        std::cout << "Order of the derivative for interpolation: " << derivativeOrder <<
//                     " not supported" << std::endl;
//        exit(EXIT_FAILURE);
//    }

//    // obtain derivative approximation m
//    for (int i=0; i<dimension; i++){
//        pk = n;
//        pkp1 = n; pkp1[i] += 1;
//        pkm1 = n; pkm1[i] -= 1;
//        vector<double> mArray;
//        mArray.push_back(0.5*((getEval(pkp1) - getEval(pk))/
//                       (discretization[i][n[i]+1] - discretization[i][n[i]]) +
//                       (getEval(pk) - getEval(pkm1))/
//                       (discretization[i][n[i]] - discretization[i][n[i]-1])));

//        pkp2 = pkp1; pkp2[i] += 1;
//        mArray.push_back(0.5*((getEval(pkp2) - getEval(pkp1))/
//                       (discretization[i][n[i]+2] - discretization[i][n[i]+1]) +
//                       (getEval(pkp1) - getEval(n))/
//                       (discretization[i][n[i]+1] - discretization[i][n[i]])));

//        m.push_back(mArray);
//    }

//    for (int i=0; i<dimension; i++){
//        Beta0 *= beta[i][0];
//        Beta1 *= beta[i][1];
//        Beta2 *= beta[i][2];
//        Beta3 *= beta[i][3];
//    }
////    std::cout << "Beta0: " << Beta0 << std::endl;
////    std::cout << "Beta1: " << Beta1 << std::endl;
////    std::cout << "Beta2: " << Beta2 << std::endl;
////    std::cout << "Beta3: " << Beta3 << std::endl;

//    for (int i=0; i<dimension; i++){
//        pkp1 = n; pkp1[i] += 1;
//        z += (Beta0*getEval(n) + Beta1*m[i][0] +
//              Beta2*getEval(pkp1) + Beta3*m[i][1]);
////        z += (beta[i][0]*getEval(n) + beta[i][1]*m[i][0] +
////              beta[i][2]*getEval(pkp1) + beta[i][3]*m[i][1]);
////        std::cout << "\n" << "z[" << i << "]: " <<
////                     (Beta0*getEval(n) + Beta1*m[i][0] + Beta2*getEval(pkp1) + Beta3*m[i][1]) << std::endl;
//    }
//    return z;
//}



//void Interpolate::computeBasisFunctionsDerivatives(vector<vector<double>> &beta,
//                                              vector<double> u,
//                                              int order){
////    // compute the derivatives of the basis spline function for evaluation of the
////    // derivative of the grid-frield
////    for (int i=0; i<dimension; i++){
////        double Bk[4];
////        for (int k=0; k<4; k++){
////            Bk[k] = binomialCoefficient(order,k)*k*pow(u[i],k-1)*(order-k)*pow(1-u[i],order-k-1);
////        }
////        vector<double> betaArray;
////        betaArray.push_back(Bk[0] + Bk[1]);
////        betaArray.push_back(Bk[1]/3);
////        betaArray.push_back(Bk[3] + Bk[2]);
////        betaArray.push_back(-Bk[2]/3);
////        beta.push_back(betaArray);
////    }
//    for (int i=0; i<dimension; i++){
//        vector<double> betaArray;
//        betaArray.push_back(6*pow(u[i],2) - 6*u[i]);
//        betaArray.push_back(3*pow(u[i],2) - 4*u[i] + 1);
//        betaArray.push_back(-6*pow(u[i],2) + 6*u[i]);
//        betaArray.push_back(3*pow(u[i],2) - 2*u[i]);
//        beta.push_back(betaArray);
//    }
//}



//    for (int i=0; i<dimension; i++){
//        vector<double> betaArray;
//        betaArray.push_back(2*pow(u[i],3) - 3*pow(u[i],2) + 1);
//        betaArray.push_back(pow(u[i],3) - 2*pow(u[i],2) + u[i]);
//        betaArray.push_back(-2*pow(u[i],3) + 3*pow(u[i],2));
//        betaArray.push_back(pow(u[i],3) - pow(u[i],2));
//        beta.push_back(betaArray);
//    }

//    for (int i=0; i<dimension; i++){
//        vector<double> betaArray;
//        betaArray.push_back(6*pow(u[i],2) - 6*u[i]);
//        betaArray.push_back(3*pow(u[i],2) - 4*u[i] + 1);
//        betaArray.push_back(-6*pow(u[i],2) + 6*u[i]);
//        betaArray.push_back(3*pow(u[i],2) - 2*u[i]);
//        beta.push_back(betaArray);
//    }








//    // get the index of the closest range value to the discretization point
//    for (int i=0; i<dimension; i++){
//        auto it = std::find_if(std::begin(discretization[i]),
//                               std::end(discretization[i]),
//                               [&](double j){return j > x[i];});
//        n[i] = std::distance(discretization[i].begin(), it)-1;
//    }

//    // compute remaining fraction
//    for (int i=0; i<dimension; i++){
//        u[i] = (x[i]-discretization[i][n[i]])/
//                (discretization[i][2]-discretization[i][1]);
//    }

//    // compute the polynomials (already evaluated)
//    for (int i=0; i<dimension; i++){
//        // compute binomial coefficient derivatives
//        beta[i][0] = 5*pow(u[i],4) - 10*pow(u[i],3) + 4.5*pow(u[i],2) + u[i] - 0.5;
//        beta[i][1] = -15*pow(u[i],4) + 30*pow(u[i],3) - 13.5*pow(u[i],2) - 2*u[i];
//        beta[i][2] = 15*pow(u[i],4) - 30*pow(u[i],3) + 13.5*pow(u[i],2) + u[i] + 0.5;
//        beta[i][3] = pow(u[i],2)*(-5*pow(u[i],2) + 10*u[i] - 4.5);
//    }

//    // loop over all the considered points (n-dimensional) and multiply the
//    // evaluation with the weight
//    double z, Beta;
//    int discrLoopCnt = -1;
//    for (int i=0; i<4; i++){
//        Beta = beta[coordinate][i];
//        for (int j=0; j<dimension; j++){
//            loc[j] = n[j];
//        }
//        loc[coordinate] += discrLoopCnt;

//        // this division (/) works for 1D but for nD things need to change
//        z += getEvalFast()*Beta/
//                (discretization[coordinate][1]-discretization[coordinate][0]);

//        discrLoopCnt += 1;
//    }



//    int discrLoopCnt[dimension] = {-1};
//    double z;
//    bool breakWhile = false;
//    bool allTrue;
//    double Beta = 1;

//    while (discrLoopCnt[0] < 3){
//        Beta = 1;
//        for (int i=0; i<dimension; i++){
//            Beta = Beta*beta[i][discrLoopCnt[i]+g];
//        }

//        for (int i=0; i<dimension; i++){
//            loc[i] = discrLoopCnt[i] + n[i];
//        }
//        z += getEvalFast()*Beta;


//        // EXIT CONDITIONS
//        for (int x=dimension-1; x>=0; x--){
//            if (discrLoopCnt[x] != 2){
//                discrLoopCnt[x] += 1;
//                break;
//            }
//            if (discrLoopCnt[x] == 2){
//                for (int y=x; y<dimension; y++){
//                    discrLoopCnt[y] = -1;
//                }
//            }
//        }
//        if (breakWhile){
//            break;
//        }
//        allTrue = true;
//        for (int i=0; i<dimension; i++){
//            if (discrLoopCnt[i] != 2){
//                allTrue = false;
//                break;
//            }
//        }
//        if (allTrue){
//            breakWhile = true;
//        }
//    }
//    return z;
