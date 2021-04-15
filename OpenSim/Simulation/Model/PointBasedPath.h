#ifndef OPENSIM_POINTBASED_PATH_H_
#define OPENSIM_POINTBASED_PATH_H_

/* -------------------------------------------------------------------------- *
 *                          OpenSim:  GeometryPath.h                          *
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

#ifdef SWIG
    #ifdef OSIMSIMULATION_API
        #undef OSIMSIMULATION_API
        #define OSIMSIMULATION_API
    #endif
#endif

// Forward declarations
namespace Simbody {
class State;
}

namespace OpenSim {
class Coordinate;
}

namespace OpenSim {

/**
 * An `OpenSim::GeometryPath` that is implemented as a sequence of points
 * connected to their adjacent neighbors.
 *
 * Historically, in earlier versions of OpenSim, `OpenSim::GeometryPath`s were
 * always implemented as `PointBasedPath`s. In later versions, there was a
 * desire to define paths in terms of mathematical functions, curve fitting,
 * lookups, etc., which is why the distinction now exists.
 */
class OSIMSIMULATION_API PointBasedPath : public GeometryPath {
    OpenSim_DECLARE_CONCRETE_OBJECT(PointBasedPath, GeometryPath);

public:
    PointBasedPath();

    double getLength(const SimTK::State& s) const override;
    void setLength(const SimTK::State& s, double length) const override;

    double getLengtheningSpeed(const SimTK::State& s) const override;
    void setLengtheningSpeed(const SimTK::State& s, double speed) const override;

    double computeMomentArm(const SimTK::State& s,
                            const Coordinate& aCoord) const override;
};

}
#endif
