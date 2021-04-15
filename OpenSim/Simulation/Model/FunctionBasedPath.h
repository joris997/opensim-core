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
    static FunctionBasedPath fromDataFile(const std::string& path);

    FunctionBasedPath();
    FunctionBasedPath(const FunctionBasedPath&);
    FunctionBasedPath(FunctionBasedPath&&);
    FunctionBasedPath& operator=(const FunctionBasedPath&);
    FunctionBasedPath& operator=(FunctionBasedPath&&);
    ~FunctionBasedPath() noexcept override;

    double getLength(const SimTK::State& s) const override;
    void setLength(const SimTK::State& s, double length) const override;

    double getLengtheningSpeed(const SimTK::State& s) const override;
    void setLengtheningSpeed(const SimTK::State& s, double speed) const override;

    const std::string& getDataPath() const { return get_data_path(); }
    void setDataPath(const std::string& newPath) { upd_data_path() = newPath; }

    double computeMomentArm(const SimTK::State& s,
                            const Coordinate& aCoord) const override;

    void printContent(std::ostream& printFile) const;
};
}

#endif
