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

#include "PointBasedPath.h"
#include <OpenSim/Simulation/SimbodyEngine/Coordinate.h>

////////////////////
// FUNCTIONS ADAM //
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
// factory function to Defer
template<typename Callback>
Defer<Callback> defer_action(Callback cb);

bool coord_affects_muscle(
    OpenSim::PointBasedPath const& pbp,
    OpenSim::Coordinate const& c,
    SimTK::State& state,
    Nonzero_conditions& out);

//////////////////////////
// ADDITIONAL FUNCTIONS //
//////////////////////////
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

class Interpolate{
    private:
        // dimension of the interpolation (e.g. two states -> 2D)
        int dimension;
        // number of points considered in constructing the polynomial
        int nInterPoints = 4;
        int g = 1;
        // vector containing the orthogonal discretization vectors
        // e.g. disc[0] = xrange, disc[1] = yrange etc.
        std::vector<std::vector<double>> discretization;
        std::vector<double> discSizes;
        // array containing the evaluations over the discretizations
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
        std::vector<double> getRange(int i) const {return discretization[i];}
        int getDimension() const {return dimension;}
        std::vector<Discretization> getdS() const {return dS;}
//        std::vector<OpenSim::Coordinate *> getCoords() const {return coords;}
        std::vector<double> getEvals() const {return evals;}

        // template member functions
        double getEval();

        double getInterp(const std::vector<double>& x);
        double getInterp(const SimTK::State& s);

        double getInterpStruct(const std::vector<double>& x);
        double getInterpStruct(const SimTK::State& s);
        double getLength(const SimTK::State& s){return getInterpStruct(s);}

        double getInterpDer(const std::vector<double>& x, int coordinate, double h=0.0001);
        double getInterpDer(const SimTK::State& s, int coordinate, double h=0.0001);
        double getLengtheningSpeed(const SimTK::State& s);

        // EMPTY CONSTRUCTOR
        explicit Interpolate();

        // READ DATA CONSTRUCTOR
        explicit Interpolate(std::vector<OpenSim::Coordinate const*> coordsIn,
                             std::vector<Discretization> dSIn,
                             std::vector<double> evalsIn);

        // CONSTRUCTOR FOR A NON-GEOMETRYPATH/COORDINATE INTERPOLATION
        explicit Interpolate(std::vector<std::vector<double>> discretizationIn,
                             std::vector<std::pair<std::vector<int>,double>> evalsPair);

        // N COORDINATES    WITH GIVEN NPOINTS
        explicit Interpolate(OpenSim::PointBasedPath const& pbp,
                             OpenSim::Coordinate const** cBegin,
                             OpenSim::Coordinate const** cEnd,
                             SimTK::State& st,
                             std::vector<int>& discretizationNPoints);


        // VECTOR OF COORDINATES
//        explicit Interpolate(OpenSim::Muscle const& m,
//                             std::vector<OpenSim::Coordinate const*> coords,
//                             SimTK::State& st,
//                             std::vector<int>& discretizationNPoints);

        explicit Interpolate(OpenSim::PointBasedPath const& pbp,
                             std::vector<OpenSim::Coordinate const*> coords,
                             SimTK::State& st,
                             std::vector<int>& discretizationNPoints);

};

#endif // INTERPOLATE_H
