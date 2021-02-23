#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <functional>

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimulationUtilities.h>

template<typename T>
void linspace(std::vector<double> &linspaced, T start_in, T end_in, int num_in){
// from: https://stackoverflow.com/questions/27028226/python-linspace-in-c
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

class Interpolate{
    private:
        // dimension of the interpolation (e.g. two states -> 2D)
        int dimension;
        // number of points considered in constructing the polynomial
        int nInterPoints = 4;
        int g = 1;
        // order of the interpolation polynomial function
//        int polyOrder = 5;
        // vector containing the orthogonal discretization vectors
        // e.g. disc[0] = xrange, disc[1] = yrange etc.
        std::vector<std::vector<double>> discretization;
        std::vector<double> discSizes;
        // array containing the evaluations over the discretizations
        std::vector<std::pair<std::vector<int>,double>> evalsPair;
        std::vector<double> evals;
        // vector containing index closest to discretization point
        std::vector<int> n;
        // vector containing fraction of closest index to discretization point to next
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
        std::vector<double> getRange(int i) {return discretization[i];}
        int getDimension() {return dimension;}

        // template member functions
        double getEval();
        double getEvalFast();

        double getInterp(const std::vector<double>& x);
        double getInterp(const SimTK::State& s);

        double getInterpStruct(const std::vector<double>& x);
        double getInterpStruct(const SimTK::State& s);
        double getLength(const SimTK::State& s){return getInterpStruct(s);}

        double getInterpDer(const std::vector<double>& x, int coordinate, double h=0.0001);
        double getInterpDer(const SimTK::State& s, int coordinate, double h=0.0001);
        double getLengtheningSpeed(const SimTK::State& s);

        double interpCubicHermiteSpline(std::vector<double> x, int derivativeOrder);
        void computeBasisFunctions(std::vector<double> u, int order);
        void computeBasisFunctionsDerivatives(std::vector<std::vector<double>> &beta,
                                              std::vector<double> u, int order);


        // CONSTRUCTOR FOR A NON-MUSCLE/COORDINATE INTERPOLATION
        explicit Interpolate(std::vector<std::vector<double>> discretizationIn,
                             std::vector<std::pair<std::vector<int>,double>> evalsPair);

        // N COORDINATES    WITH GIVEN NPOINTS
        explicit Interpolate(OpenSim::Muscle const& m,
                             OpenSim::Coordinate const** cBegin,
                             OpenSim::Coordinate const** cEnd,
                             SimTK::State& st,
                             std::vector<int>& discretizationNPoints);


        // VECTOR OF COORDINATES
        explicit Interpolate(OpenSim::Muscle const& m,
                             std::vector<OpenSim::Coordinate const*> coords,
                             SimTK::State& st,
                             std::vector<int>& discretizationNPoints);

//        // 1 COORDINATE
//        explicit Interpolate(OpenSim::Muscle const& m,
//                        Coordinate const& c1,
//                        SimTK::State& st,
//                        vector<int>& discretizationNPoints,
//                        bool smooth=false)
//            : Interpolate(m,
//                     c1,
//                     c1,
//                     st,
//                     discretizationNPoints) {}
};

#endif // INTERPOLATE_H







//        // 1 COORDINATE
//        explicit Interpolate(OpenSim::Muscle const& m,
//                        Coordinate const& c1,
//                        SimTK::State& st,
//                        vector<int>& discretizationNPoints,
//                        bool smooth=false)
//            : Interpolate(m,
//                     c1,
//                     c1,
//                     st,
//                     discretizationNPoints) {}

//        // 2 COORDINATES
//        explicit Interpolate(OpenSim::Muscle const& m,
//                        Coordinate const& c1,
//                        Coordinate const& c2,
//                        SimTK::State& st,
//                        vector<int>& discretizationNPoints,
//                        bool smooth=false)
//            : Interpolate(m,
//                     c1,
//                     c2,
//                     st,
//                     discretizationNPoints) {}

//        // 3 COORDINATES
//        explicit Interpolate(OpenSim::Muscle const& m,
//                        Coordinate const& c1,
//                        Coordinate const& c2,
//                        Coordinate const& c3,
//                        SimTK::State& st,
//                        vector<int>& discretizationNPoints,
//                        bool smooth=false)
//            : Interpolate(m,
//                     c1,
//                     c3,
//                     st,
//                     discretizationNPoints) {}

//        // 4 COORDINATES
//        explicit Interpolate(OpenSim::Muscle const& m,
//                        Coordinate const& c1,
//                        Coordinate const& c2,
//                        Coordinate const& c3,
//                        Coordinate const& c4,
//                        SimTK::State& st,
//                        vector<int>& discretizationNPoints,
//                        bool smooth=false)
//            : Interpolate(m,
//                     c1,
//                     c4,
//                     st,
//                     discretizationNPoints) {}




//        // N COORDINATES    WITH GIVEN DISCRETIZATION STRUCT
//        explicit Interpolate(OpenSim::Muscle const& m,
//                        Coordinate const** cBegin,
//                        Coordinate const** cEnd,
//                        SimTK::State& st,
//                        vector<Discretization>& discretizationData)
//            : dimension(discretizationData.size()),
//              n(dimension,0),
//              u(dimension,0),
//              loc(dimension,0),
//              dS(discretizationData)
//        {
//            std::ptrdiff_t n = cEnd-cBegin;
//            assert(discretizationData.size() == (int)n);

//            GeometryPath const& musc_path = m.getGeometryPath();

//            // just make it for using the old getInterp method
//            vector<double> dc_;
//            for (int i=0; i<dimension; i++){
//                beta.push_back({0,0,0,0});

//                dc_.clear();
//                linspace(dc_,dS[i].begin,dS[i].end,dS[i].nPoints);
//                discretization.push_back(dc_);
//                discSizes.push_back(dS[i].nPoints);
//            }

//            // unlock coordinates
//            for (int i=0; i<dimension; i++){
//                Coordinate c = **(cBegin+i);
//                bool c_was_locked = c.getLocked(st);
//                c.setLocked(st, false);
//                auto unlock_c = defer_action([&] { c.setLocked(st, c_was_locked); });
//                double c_initial_value = c.getValue(st);
//                auto reset_c_val = defer_action([&] { c.setValue(st, c_initial_value); });
//            }

//            int numOfLoops = 1;
//            for (int i=0; i<dimension; i++){
//                numOfLoops *= dS[i].nPoints;
//            }
//            vector<int> cnt[dimension];
//            vector<int> coordValues[dimension];
//            for (int i=0; i<numOfLoops; i++){
//                for (int ii=0; ii<dimension; ii++){
//                    coordValues[ii] = dS[ii].begin + (cnt[ii]*dS[ii].gridsize);
//                    **(cBegin+ii).setValue(st,coordValues[ii]);
//                }
//                evals.push_back((musc_path.getLength(st)));
//                // update cnt values
//                for (int x=dimension-1; x>=0; x--){
//                    if (cnt[x] != dS[x].nPoints-1){
//                        cnt[x] += 1;
//                        break;
//                    }
//                    if (cnt[x] == dS[x].nPoints-1){
//                        for (int y=x; y<dimension; y++){
//                            cnt[y] = 0;
//                        }
//                    }
//                }

//            }
//        }





















//// 1 COORDINATE
//explicit Interpolate(OpenSim::Muscle const& m,
//                Coordinate const& c1,
//                SimTK::State& st,
//                vector<Discretization>& discretizationData)
//    : dimension(discretizationData.size()),
//      n(dimension,0),
//      u(dimension,0),
//      loc(dimension,0),
//      dS(discretizationData)
//{
//    assert(discretizationData.size() == 1);
//    coords.push_back(c1);
//    GeometryPath const& musc_path = m.getGeometryPath();

//    // just make it for using the old getInterp method
//    vector<double> dc_;
//    for (int i=0; i<dimension; i++){
//        beta.push_back({0,0,0,0});

//        dc_.clear();
//        linspace(dc_,dS[i].begin,dS[i].end,dS[i].nPoints);

//        assert(dc_[1]-dc_[0] == dS[i].gridsize);

//        discretization.push_back(dc_);
//        discSizes.push_back(dS[i].nPoints);
//    }

//    // unlock coordinates
//    bool c1_was_locked = c1.getLocked(st);
//    c1.setLocked(st, false);
//    auto unlock_c1 = defer_action([&] { c1.setLocked(st, c1_was_locked); });
//    double c1_initial_value = c1.getValue(st);
//    auto reset_c1_val = defer_action([&] { c1.setValue(st, c1_initial_value); });

//    // evalute the muscle length
//    Discretization dc1 = dS[0];
//    for (int ic1=0; ic1<dc1.nPoints; ++ic1){
//        double c1v = dc1.begin + (ic1*dc1.gridsize);
//        c1.setValue(st, c1v);
//        evals.push_back(musc_path.getLength(st));
//    }
//}

//// 2 COORDINATES
//explicit Interpolate(OpenSim::Muscle const& m,
//                Coordinate const& c1,
//                Coordinate const& c2,
//                SimTK::State& st,
//                vector<Discretization>& discretizationData)
//    : dimension(discretizationData.size()),
//      n(dimension,0),
//      u(dimension,0),
//      loc(dimension,0),
//      dS(discretizationData)
//{
//    assert(discretizationData.size() == 2);
//    coords.push_back(c1);
//    coords.push_back(c2);
//    GeometryPath const& musc_path = m.getGeometryPath();

//    // just make it for using the old getInterp method
//    vector<double> dc_;
//    for (int i=0; i<dimension; i++){
//        beta.push_back({0,0,0,0});

//        dc_.clear();
//        linspace(dc_,dS[i].begin,dS[i].end,dS[i].nPoints);

//        assert(dc_[1]-dc_[0] == dS[i].gridsize);

//        discretization.push_back(dc_);
//        discSizes.push_back(dS[i].nPoints);
//    }

//    // unlock coordinates
//    bool c1_was_locked = c1.getLocked(st);
//    c1.setLocked(st, false);
//    auto unlock_c1 = defer_action([&] { c1.setLocked(st, c1_was_locked); });
//    double c1_initial_value = c1.getValue(st);
//    auto reset_c1_val = defer_action([&] { c1.setValue(st, c1_initial_value); });

//    bool c2_was_locked = c2.getLocked(st);
//    c2.setLocked(st, false);
//    auto unlock_c2 = defer_action([&] { c2.setLocked(st, c2_was_locked); });
//    double c2_initial_value = c2.getValue(st);
//    auto reset_c2_val = defer_action([&] { c2.setValue(st, c2_initial_value); });
//    // evalute the muscle length

//    Discretization dc1 = dS[0];
//    Discretization dc2 = dS[1];
//    for (int ic1=0; ic1<dc1.nPoints; ++ic1){
//        double c1v = dc1.begin + (ic1*dc1.gridsize);
//        c1.setValue(st, c1v);
//        for (int ic2=0; ic2<dc2.nPoints; ++ic2){
//            double c2v = dc2.begin + (ic2*dc2.gridsize);
//            c2.setValue(st, c2v);
//            evals.push_back(musc_path.getLength(st));
//        }
//    }
//}

//// 2 COORDINATES GIVEN N DISC POINTS
//explicit Interpolate(OpenSim::Muscle const& m,
//                Coordinate const& c1,
//                Coordinate const& c2,
//                SimTK::State& st,
//                vector<int>& discretizationNPoints,
//                bool smooth=false)
//    : dimension(discretizationNPoints.size()),
//      n(dimension,0),
//      u(dimension,0),
//      loc(dimension,0)
//{
//    assert(discretizationNPoints.size() == 2);
//    coords.push_back(c1);
//    coords.push_back(c2);
//    std::cout << "-------- CONSTRUCTOR CALLED --------" << std::endl;
//    GeometryPath const& musc_path = m.getGeometryPath();

//    // unlock coordinates
//    bool c1_was_locked = c1.getLocked(st);
//    c1.setLocked(st, false);
//    auto unlock_c1 = defer_action([&] { c1.setLocked(st, c1_was_locked); });
//    double c1_initial_value = c1.getValue(st);
//    auto reset_c1_val = defer_action([&] { c1.setValue(st, c1_initial_value); });

//    bool c2_was_locked = c2.getLocked(st);
//    c2.setLocked(st, false);
//    auto unlock_c2 = defer_action([&] { c2.setLocked(st, c2_was_locked); });
//    double c2_initial_value = c2.getValue(st);
//    auto reset_c2_val = defer_action([&] { c2.setValue(st, c2_initial_value); });

//    // make discretization objects for interpolation class instance
//    Discretization dc1, dc2;
//    dc1.begin    = c1.getRangeMin();
//    dc1.end      = c1.getRangeMax();
//    dc1.nPoints  = discretizationNPoints[0];
//    dc1.gridsize = (dc1.end-dc1.begin) / (dc1.nPoints-1);
//    dc2.begin    = c2.getRangeMin();
//    dc2.end      = c2.getRangeMax();
//    dc2.nPoints  = discretizationNPoints[1];
//    dc2.gridsize = (dc2.end-dc2.begin) / (dc2.nPoints-1);
//    // push into discretization property vector
//    dS.push_back(dc1); dS.push_back(dc2);

//    // check rms of error and possibly refine grid
//    if (smooth){
//        std::cout << "-------- STARTING GRID SHAPING --------" << std::endl;
//        // make it 6 points so we have enough room to sample inside. If it's
//        // too large from user's input we throw away computation time, if it's
//        // too small we cannot create sample points of which the surrounded
//        // area does not exceed the min max of the coordinate
//        for (int dim=0; dim<dimension; dim++){
//            dS[dim].nPoints = 10;
//            dS[dim].gridsize = (dS[dim].end-dS[dim].begin)/(dS[dim].nPoints-1);
//        }
//        double RMS;
//        int whileCnt;
//        // this section samples 'num_of_points' number of points on the existing
//        // grid and creates an interpolation object with bounds around it that
//        // never exceed the min and max of the coordinate. It then samples a point
//        // inside this interpolation object and evaluates the value. This is compared
//        // to the real value and the gridsize is refined if the error remains too large
//        while (whileCnt < 6){
//            // create a grid on which we would like to sample for accuracy
//            vector<vector<double>> random_coords;
//            vector<double> samples;
//            int num_of_points = 10;
//            for (int sample=0; sample<num_of_points+1; sample++){
//                // create a vector of samples for a certain coordinate
//                samples.clear();
//                for (int dim=0; dim<dimension; dim++){
//                    std::uniform_real_distribution<double> unif(
//                                dS[dim].begin+dS[dim].gridsize,
//                                dS[dim].end-2*dS[dim].gridsize);
//                    std::random_device rd;
//                    std::mt19937 gen(rd());
//                    samples.push_back(unif(gen));
//                }
//                random_coords.push_back(samples);
//            }
//            // include the total edge case as 'random' point
//            samples.clear();
//            for (int dim=0; dim<dimension; dim++){
//                samples.push_back(dS[dim].end-2*dS[dim].gridsize);
//            }
//            random_coords.push_back(samples);

//            // create a mini interp object for a single point
//            vector<double> errors;
//            for (int i=0; i<random_coords.size(); i++){
//                vector<Discretization> dcVec;
//                Discretization dc;
//                for (int dim=0; dim<dimension; dim++){
//                    // consider 1 point before and 2 points after for accuracy on edges
//                    dc.begin    = random_coords[i][dim] - 1*dS[dim].gridsize;
//                    dc.end      = random_coords[i][dim] + 2*dS[dim].gridsize;
//                    dc.gridsize = dS[dim].gridsize;
//                    dc.nPoints  = 4;
//                    dcVec.push_back(dc);

//                    // get some points around the center point and evaluate
//                    std::uniform_real_distribution<double> unif(0,dc.gridsize);
//                    std::random_device rd;
//                    std::mt19937 gen(rd());
//                    // we can overwrite [i][dim] and add a sampling from begin to end
//                    // as this lays in front or behind of random_coords[i][dim]
//                    random_coords[i][dim] += unif(gen);
//                }
//                // get approx
//                Interpolate a = Interpolate(m,c1,c2,st,dcVec);
//                double lengthInterp = a.getInterpStruct(random_coords[i]);
//                std::cout << "interp: " << lengthInterp << std::endl;
//                // get real
//                c1.setValue(st,random_coords[i][0]);
//                c2.setValue(st,random_coords[i][1]);
//                double lengthReal = musc_path.getLength(st);
//                std::cout << "real:   " << lengthReal << std::endl;
//                // get error value
//                errors.push_back(std::abs(lengthInterp-lengthReal)/lengthReal);
//            }
//            RMS = *max_element(errors.begin(),errors.end());
//            std::cout << "nPoints: " << dS[0].nPoints << std::endl;
//            std::cout << "RMS: " << RMS << std::endl;
//            if (RMS < 0.001){
//                break;
//            } else {
//                for (int dim=0; dim<dimension; dim++){
//                    // double number of interpolation points
//                    dS[dim].nPoints  = 2*dS[dim].nPoints - 1;
//                    dS[dim].gridsize = (dS[dim].end-dS[dim].begin)/(dS[dim].nPoints-1);
//                }
//            }
//            whileCnt++;
//        }
//    }

//    // just make it for using the old getInterp method
//    vector<double> dc_;
//    for (int i=0; i<dimension; i++){
//        beta.push_back({0,0,0,0});

//        dc_.clear();
//        linspace(dc_,dS[i].begin,dS[i].end,dS[i].nPoints);
//        discretization.push_back(dc_);
//        discSizes.push_back(dS[i].nPoints);
//    }

//    // evalute the muscle length
//    std::cout << "-------- EVALUATING REAL LENGTHS --------" << std::endl;
//    for (int ic1=0; ic1<dc1.nPoints; ++ic1){
//        double c1v = dc1.begin + (ic1*dc1.gridsize);
//        c1.setValue(st, c1v);
//        for (int ic2=0; ic2<dc2.nPoints; ++ic2){
//            double c2v = dc2.begin + (ic2*dc2.gridsize);
//            c2.setValue(st, c2v);
//            evals.push_back(musc_path.getLength(st));
//        }
//    }
//    std::cout << "-------- CONSTRUCTOR FINISHED --------\n" << std::endl;
//}












