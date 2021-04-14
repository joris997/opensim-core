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
        explicit Interpolate();

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

#endif // INTERPOLATE_H
