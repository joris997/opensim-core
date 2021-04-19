/* -------------------------------------------------------------------------- *
 *                         OpenSim:  testForward.cpp                          *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
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
// functionbasedpathconversion.cpp
// author:  Joris Verhagen

// INCLUDE
#include <string>
#include <iostream>
#include <OpenSim/OpenSim.h>

using namespace OpenSim;
using namespace std;

void testArm();


int main() {
    SimTK_START_TEST("testFunctionBasedPathConversion");
        // test arm model conversion accuracy
        SimTK_SUBTEST(testArm);
    SimTK_END_TEST();
}

void testArm(){
    try {
        string modelName = "arm26.osim";
        string newModelName = "arm26_FBP.osim";

        // test creating a FBP model
        FunctionBasedPathConversionTool tool(modelName,newModelName);
        tool.run();

        // test loading in a FBP model
        Model modelPoint(modelName);
        Model modelFunction(newModelName);
        // load/overwrite interp from file

        modelPoint.finalizeConnections();
        modelPoint.finalizeFromProperties();
        modelPoint.initSystem();
        modelPoint.printSubcomponentInfo();

        modelFunction.finalizeConnections();
        modelFunction.finalizeFromProperties();
        modelFunction.initSystem();
        modelFunction.printSubcomponentInfo();

        // test simulating a FBP model
        auto reporterPoint = new ConsoleReporter();
        reporterPoint->setName("point_results");
        reporterPoint->set_report_time_interval(0.05);
        reporterPoint->addToReport(modelPoint.getComponent("/forceset/TRIlong/pointbasedpath").getOutput("length"));
        modelPoint.addComponent(reporterPoint);

        auto reporterFunction = new ConsoleReporter();
        reporterFunction ->setName("function_results");
        reporterFunction ->set_report_time_interval(0.05);
        reporterFunction ->addToReport(modelFunction.getComponent("/forceset/TRIlong/functionbasedpath").getOutput("length"));
        modelFunction.addComponent(reporterFunction);

        SimTK::State& pointState = modelPoint.initSystem();
        SimTK::State& functionState = modelFunction.initSystem();

        chrono::milliseconds pbpSimTime{0};
        chrono::milliseconds fbpSimTime{0};
        int pbpStepsAttempted = 0;
        int fbpStepsAttempted = 0;
        int pbpStepsTaken = 0;
        int fbpStepsTaken= 0;

        int n = 10;
        double tFinal = 3.5;
        for (int i=0; i<n; i++){
            OpenSim::Manager manager(modelPoint);
            manager.initialize(pointState);
            auto before = chrono::high_resolution_clock::now();
            manager.integrate(tFinal);
            auto after = chrono::high_resolution_clock::now();
            auto dt = after - before;
            pbpStepsAttempted += manager.getIntegrator().getNumStepsAttempted();
            pbpStepsTaken += manager.getIntegrator().getNumStepsTaken();
            pbpSimTime += chrono::duration_cast<chrono::milliseconds>(dt);
        }

        for (int i=0; i<n; i++){
            OpenSim::Manager manager(modelFunction);
            manager.initialize(functionState);
            auto before = chrono::high_resolution_clock::now();
            manager.integrate(tFinal);
            auto after = chrono::high_resolution_clock::now();
            auto dt = after - before;
            fbpStepsAttempted += manager.getIntegrator().getNumStepsAttempted();
            fbpStepsTaken += manager.getIntegrator().getNumStepsTaken();
            fbpSimTime += chrono::duration_cast<chrono::milliseconds>(dt);
        }
        pbpSimTime /= n;
        fbpSimTime /= n;
        pbpStepsAttempted /= n;
        fbpStepsAttempted /= n;
        pbpStepsTaken /= n;
        fbpStepsTaken /= n;

        cout << "\navg-time PBP:        " << pbpSimTime.count() << endl;
        cout << "avg-time FBP:        "   << fbpSimTime.count() << endl;

        cout << "\nsteps attempted PBP: " << pbpStepsAttempted << endl;
        cout << "steps attempted FBP: "   << fbpStepsAttempted << endl;

        cout << "\nsteps taken PBP:     " << pbpStepsTaken << endl;
        cout << "steps taken FBP:     "   << fbpStepsTaken << endl;

        // test accuracy of FBP approximation

    } catch (const Exception& e) {
        e.print(cerr);
    }
    cout << "Done" << endl;
}
