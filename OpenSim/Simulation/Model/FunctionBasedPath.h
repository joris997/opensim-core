#ifndef OPENSIM_FUNCTIONBASED_PATH_H_
#define OPENSIM_FUNCTIONBASED_PATH_H_

/* -------------------------------------------------------------------------- *
 *                          OpenSim: FunctionBasedPath.h                      *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2021 TU Delft and the Authors                           *
 * Author(s): Joris Verhagen                                                  *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include <OpenSim/Simulation/Model/GeometryPath.h>

#include <iosfwd>
#include <memory>
#include <string>
#include <utility>

#ifdef SWIG
    #ifdef OSIMSIMULATION_API
        #undef OSIMSIMULATION_API
        #define OSIMSIMULATION_API
    #endif
#endif

// forward declarations
namespace OpenSim {
class Model;
class PointBasedPath;
class Coordinate;
}

namespace SimTK {
class State;
}

namespace OpenSim {

/**
 * An `OpenSim::GeometryPath` that computes its length/lengtheningSpeed by
 * interpolating a pre-computed Bezier curve.
 *
 * The precomputed curve can be produced from by fitting an existing
 * `OpenSim::PointBasedPath` in a Model. You can use
 * `FunctionBasedPath::fromPointBasedPath` to do that in code, or use the
 * FunctionBasedPathModelTool do convert all `PointBasedPaths` in an osim into
 * `FunctionBasedPath`s.
 *
 * The curve fitting implementation requires access to the whole OpenSim model
 * that the (concrete) `PointBasedPath` is in. This is because the fitting
 * implementation permutes the model's coordinates to fit the curve, and
 * because a curve's paramaterization depends on its place in the Model as
 * a whole.
 */
class OSIMSIMULATION_API FunctionBasedPath : public GeometryPath {    
    OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPath, GeometryPath);

    OpenSim_DECLARE_PROPERTY(data_path, std::string, "Path to the associated data file for this path. Typically, the file contains interpolated data that the implementation uses to compute the path's length and lengthening speed at simulation-time");

    struct Impl;
private:
    SimTK::ClonePtr<Impl> _impl;

public:
    /**
     * Returns an in-memory representation of a FunctionBasedPath generated
     * by fitting curves against the supplied PointBasedPath.
     *
     * Does not set the 'data' member (it's an in-memory represenation). Saving
     * the resulting FunctionBasedPath without setting 'data' will throw an
     * exception.
     */
    static FunctionBasedPath fromPointBasedPath(const Model& model,
                                                const PointBasedPath& pbp);

    /**
     * Returns an FunctionBasedPath read from an associated data file
     *
     * The file must contain the content written by `printContent`
     */
    FunctionBasedPath fromDataFile(const std::string& path);

    FunctionBasedPath();
    FunctionBasedPath(const FunctionBasedPath&);
    FunctionBasedPath(FunctionBasedPath&&);
    FunctionBasedPath& operator=(const FunctionBasedPath&);
    FunctionBasedPath& operator=(FunctionBasedPath&&);
    ~FunctionBasedPath() noexcept override;

    void extendFinalizeFromProperties() override;
    void extendFinalizeConnections(Component &root) override;

    double getLength(const SimTK::State& s) const override;
    void setLength(const SimTK::State& s, double length) const override;

    double getLengtheningSpeed(const SimTK::State& s) const override;
    void setLengtheningSpeed(const SimTK::State& s, double speed) const override;

    const std::string& getDataPath() const { return get_data_path(); }
    void setDataPath(const std::string& newPath) { upd_data_path() = newPath; }

    double computeMomentArm(const SimTK::State& s,
                            const Coordinate& aCoord) const override;

    void printContent(std::ostream& printFile) const;


    // From GeometryPath refactoring
public:
    void addInEquivalentForces(const SimTK::State& state,
                               const double& tension,
                               SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                               SimTK::Vector& mobilityForces) const;

    void getPointForceDirections(const SimTK::State& s,
            OpenSim::Array<PointForceDirection*> *rPFDs) const;

protected:
    void computePath(const SimTK::State& s ) const;
    void computeLengtheningSpeed(const SimTK::State& s) const;

    double calcLengthAfterPathComputation
       (const SimTK::State& s, const Array<AbstractPathPoint*>& currentPath) const;

    void generateDecorations(
                    bool                                        fixed,
                    const ModelDisplayHints&                    hints,
                    const SimTK::State&                         state,
                    SimTK::Array_<SimTK::DecorativeGeometry>&   appendToThis) const
                    override;





    // related to pathpoints so has to be removed in release
protected:
    double calcPathLengthChange(const SimTK::State& s, const WrapObject& wo,
                                const WrapResult& wr,
                                const Array<AbstractPathPoint*>& path) const;

public:
    void addPathWrap(WrapObject& aWrapObject);

private:
    const Array<AbstractPathPoint*>& getCurrentPath( const SimTK::State& s) const;
    AbstractPathPoint* addPathPoint(const SimTK::State& s,
                                    int index,
                                    const PhysicalFrame& frame);
    AbstractPathPoint* appendNewPathPoint(const std::string& proposedName,
                                          const PhysicalFrame& frame,
                                          const SimTK::Vec3& locationOnFrame);
    bool canDeletePathPoint(int index);
    bool deletePathPoint(const SimTK::State& s,
                         int index);
    bool replacePathPoint(const SimTK::State& s,
                          AbstractPathPoint* oldPathPoint,
                          AbstractPathPoint* newPathPoint);

    void moveUpPathWrap(const SimTK::State& s,
                        int index);
    void moveDownPathWrap(const SimTK::State& s,
                          int index);
    void deletePathWrap(const SimTK::State& s,
                        int index);

    void applyWrapObjects(const SimTK::State& s,
                          Array<AbstractPathPoint*>& path) const;

    void namePathPoints(int aStartingIndex);
    void placeNewPathPoint(const SimTK::State& s, SimTK::Vec3& aOffset,
                           int index, const PhysicalFrame& frame);
};
}

#endif
