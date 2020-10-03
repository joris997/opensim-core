/* -------------------------------------------------------------------------- *
 *                         OpenSim:  testWrapping.cpp                         *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Ajay Seth                                                       *
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
// INCLUDE
#include <OpenSim/OpenSim.h>
#include <OpenSim/Auxiliary/auxiliaryTestFunctions.h>

#include "simbody/internal/CableTrackerSubsystem.h"
#include "simbody/internal/CablePath.h"
#include "simbody/internal/Force_Custom.h"

#include <set>
#include <string>
#include <iostream>

using namespace OpenSim;
using namespace SimTK;
using namespace std;

class TestInfo {
public:
    TestInfo(String modelFilename, double simulationDuration) :
        filename(modelFilename), duration(simulationDuration) {};
    TestInfo(String modelFilename) :
        filename(modelFilename), duration(1) {};
    String filename;
    Real duration;
};

// Struct that contains parameters for the analytical check of the cylinder wrapping surfaces
struct testCase {
    // General parameters
    bool   SHOW_VISUALIZER   { false };
    double REPORTING_INTERVAL{ 0.2 };
    double FINAL_TIME        { 5.0 };
    // Sliding body parameters
    double BODY_SIZE         { 0.1 };
    double BODY_OFFSET       { 0.4 };
    // Wrapping body parameters
    double CYLINDER_RADIUS   { 0.08 };
    double CYLINDER_HEIGHT   { 1.00 };
};

void testWrapCylinder();
void testWrapCylinderAnalytical();
double analyticalSolution(double leftHeight, const testCase& tc);
void testWrapObjectUpdateFromXMLNode30515();
void simulate(Model& osimModel, State& si, double initialTime, double finalTime);
void simulateModelWithMusclesNoViz(const string &modelFile, double finalTime, double activation=0.5);
void simulateModelWithPassiveMuscles(const string &modelFile, double finalTime);
void simulateModelWithoutMuscles(const string &modelFile, double finalTime);
void simulateModelWithLigaments(const string &modelFile, double finalTime);
void simulateModelWithCables(const string &modelFile, double finalTime);

int main()
{
    SimTK::Array_<std::string> failures;

    try{
        testWrapCylinder();
        // performance of multiple paths with wrapping in upper-extremity
        simulateModelWithMusclesNoViz("TestShoulderWrapping.osim", 0.1);}
    catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
        failures.push_back("TestShoulderModel (multiple wrap)"); }

    try{
        testWrapCylinderAnalytical();
    } catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
        failures.push_back("test cylinder wrap with analytical solution");
    }

    try{
        testWrapObjectUpdateFromXMLNode30515();
    } catch (const std::exception& e) {
         std::cout << "Exception: " << e.what() << std::endl;
         failures.push_back("testWrapObjectUpdateFromXMLNode30515");
    }

    if (!failures.empty()) {
        cout << "Done, with failure(s): " << failures << endl;
        return 1;
    }

    cout << "Done" << endl;
    return 0;
}

int main_new()
{
    SimTK::Array_<TestInfo> tests;
    tests.push_back(TestInfo("test_wrapCylinder_vasint.osim"));
    tests.push_back(TestInfo("test_wrapEllipsoid_vasint.osim"));
    tests.push_back(TestInfo("arm26_crop.osim"));
//    tests.push_back(TestInfo("gait2392_pelvisFixed.osim",0.05));
//    tests.push_back(TestInfo("Arnold2010_pelvisFixed_crop.osim"));
//    tests.push_back(TestInfo("TestShoulderModel_crop.osim"));
//    tests.push_back(TestInfo("Arnold2010_pelvisFixed.osim"));
//    tests.push_back(TestInfo("TestShoulderModel.osim"));

    SimTK::Array_<String> failures;

    for (unsigned int i = 0; i < tests.size(); ++i) {
        cout << "testing " << tests[i].filename << " for " << tests[i].duration << " s" << endl;
        try { // performance with cylinder wrapping
            simulateModelWithPassiveMuscles(tests[i].filename, tests[i].duration);
            simulateModelWithCables(tests[i].filename, tests[i].duration);

        } catch (const std::exception& e) {
            std::cout << "Exception: " << e.what() << std::endl;
            failures.push_back(tests[i].filename);
        }
    }

    if (!failures.empty()) {
        cout << "Done, with failure(s): " << failures << endl;
        return 1;
    }

    cout << "Done" << endl;
    return 0;
}


class ObstacleInfo {
public:
    String bodyName;
    Transform X_BS; // X_BS.p used for via point location
    const WrapObject* wrapObjectPtr;
    Vec3 P_S;
    Vec3 Q_S;
    bool isVia;
    bool isActive;
};

class CableInfo {
public:
    String orgBodyName;
    Vec3 orgLoc;
    String insBodyName;
    Vec3 insLoc;

    Array<ObstacleInfo> obstacles;
};

// This force element implements an elastic cable of a given nominal length,
// and a stiffness k that generates a k*x force opposing stretch beyond
// nominal. There is also a damping term c*xdot that applies only when the
// cable is stretched and is being extended (x>0 && xdot>0). We keep track
// of dissipated power here so we can use conservation of energy to check that
// the cable and force element aren't obviously broken.
class MyCableSpringImpl : public SimTK::Force::Custom::Implementation {
public:
    MyCableSpringImpl(const GeneralForceSubsystem& forces,
                      const CablePath& path,
                      Real stiffness, Real nominal, Real damping)
    :   forces(forces), path(path), k(stiffness), x0(nominal), c(damping)
    {   ASSERT(stiffness >= 0 && nominal >= 0 && damping >= 0); }

    const CablePath& getCablePath() const {return path;}

    // Must be at stage Velocity. Evaluates tension if necessary.
    Real getTension(const State& state) const {
        ensureTensionCalculated(state);
        return Value<Real>::downcast(forces.getCacheEntry(state, tensionx));
    }

    // Must be at stage Velocity.
    Real getPowerDissipation(const State& state) const {
        const Real stretch = calcStretch(state);
        if (stretch == 0) return 0;
        const Real rate = path.getCableLengthDot(state);
        return k*stretch*std::max(c*rate, -1.)*rate;
    }

    // This integral is always available.
    Real getDissipatedEnergy(const State& state) const {
        return forces.getZ(state)[workx];
    }

    //--------------------------------------------------------------------------
    //                       Custom force virtuals

    // Ask the cable to apply body forces given the tension calculated here.
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   override
    {   path.applyBodyForces(state, getTension(state), bodyForces); }

    // Return the potential energy currently stored by the stretch of the cable.
    Real calcPotentialEnergy(const State& state) const override {
        const Real stretch = calcStretch(state);
        if (stretch == 0) return 0;
        return k*square(stretch)/2;
    }

    // Allocate the state variable for tracking dissipated energy, and a
    // cache entry to hold the calculated tension.
    void realizeTopology(State& state) const override {
        Vector initWork(1, 0.);
        workx = forces.allocateZ(state, initWork);
        tensionx = forces.allocateLazyCacheEntry(state, Stage::Velocity,
                                             new Value<Real>(NaN));
    }

    // Report power dissipation as the derivative for the work variable.
    void realizeAcceleration(const State& state) const override {
        Real& workDot = forces.updZDot(state)[workx];
        workDot = getPowerDissipation(state);
    }
    //--------------------------------------------------------------------------

private:
    // Return the amount by which the cable is stretched beyond its nominal
    // length or zero if the cable is slack. Must be at stage Position.
    Real calcStretch(const State& state) const {
        const Real stretch = path.getCableLength(state) - x0;
        return std::max(stretch, 0.);
    }

    // Must be at stage Velocity to calculate tension.
    Real calcTension(const State& state) const {
        const Real stretch = calcStretch(state);

//        cout << "cable stretch = " << stretch << endl;

        if (stretch == 0) return 0;
        const Real rate = path.getCableLengthDot(state);
        if (c*rate < -1)
            cout << "c*rate=" << c*rate << "; limited to -1\n";
        const Real tension = k*stretch*(1+std::max(c*rate,-1.));
        return tension;
    }

    // If state is at stage Velocity, we can calculate and store tension
    // in the cache if it hasn't already been calculated.
    void ensureTensionCalculated(const State& state) const {
        if (forces.isCacheValueRealized(state, tensionx))
            return;
        Value<Real>::updDowncast(forces.updCacheEntry(state, tensionx))
            = calcTension(state);
        forces.markCacheValueRealized(state, tensionx);
    }

    const GeneralForceSubsystem&    forces;
    CablePath                       path;
    Real                            k, x0, c;
    mutable ZIndex                  workx;
    mutable CacheEntryIndex         tensionx;
};

// A nice handle to hide most of the cable spring implementation. This defines
// a user's API.
class MyCableSpring : public SimTK::Force::Custom {
public:
    MyCableSpring(GeneralForceSubsystem& forces, const CablePath& path,
                  Real stiffness, Real nominal, Real damping)
    :   SimTK::Force::Custom(forces, new MyCableSpringImpl(forces,path,
                                                    stiffness,nominal,damping))
    {}

    // Expose some useful methods.
    const CablePath& getCablePath() const
    {   return getImpl().getCablePath(); }
    Real getTension(const State& state) const
    {   return getImpl().getTension(state); }
    Real getPowerDissipation(const State& state) const
    {   return getImpl().getPowerDissipation(state); }
    Real getDissipatedEnergy(const State& state) const
    {   return getImpl().getDissipatedEnergy(state); }

private:
    const MyCableSpringImpl& getImpl() const
    {   return dynamic_cast<const MyCableSpringImpl&>(getImplementation()); }
};

static Array_<State> saveStates;
// This gets called periodically to dump out interesting things about
// the cables and the system as a whole.
class ShowStuff : public PeriodicEventReporter {
public:

    ShowStuff(const MultibodySystem& mbs,
              Real interval)
    :   PeriodicEventReporter(interval) {}

    static void showHeading(std::ostream& o) {
        printf("%8s %10s %10s %10s %10s %10s %10s %10s %10s %12s\n",
            "time", "length", "rate", "integ-rate", "unitpow", "tension", "disswork",
            "KE", "PE", "KE+PE-W");
    }

    /** This is the implementation of the EventReporter virtual. **/
    void handleEvent(const State& state) const override {
//        const CablePath& path1 = cable1.getCablePath();
////        const CablePath& path2 = cable2.getCablePath();
//        printf("%8g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g CPU=%g\n",
//            state.getTime(),
//            path1.getCableLength(state),
//            path1.getCableLengthDot(state),
//            path1.getIntegratedCableLengthDot(state),
//            path1.calcCablePower(state, 1), // unit power
//            cable1.getTension(state),
//            cable1.getDissipatedEnergy(state),
//            cpuTime());
//        printf("%8s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %12.6g\n",
//            "",
//            path2.getCableLength(state),
//            path2.getCableLengthDot(state),
//            path2.getIntegratedCableLengthDot(state),
//            path2.calcCablePower(state, 1), // unit power
//            cable2.getTension(state),
//            cable2.getDissipatedEnergy(state),
//            mbs.calcKineticEnergy(state),
//            mbs.calcPotentialEnergy(state),
//            mbs.calcEnergy(state)
//                + cable1.getDissipatedEnergy(state)
//                + cable2.getDissipatedEnergy(state));
        saveStates.push_back(state);
    }
private:
//    MyCableSpring           cable1;
};


void testWrapCylinder()
{
    const double r = 0.25;
    const double off = sqrt(2)*r-0.05;
    Model model;
    model.setName("testWrapCylinder");

    auto& ground = model.updGround();
    auto body = new OpenSim::Body("body", 1, Vec3(0), Inertia(0.1, 0.1, 0.01));
    model.addComponent(body);

    auto bodyOffset = new PhysicalOffsetFrame("bToj", *body, Transform(Vec3(-off, 0, 0)));
    model.addComponent(bodyOffset);
    
    auto joint = new PinJoint("pin", ground, *bodyOffset);
    model.addComponent(joint);

    WrapCylinder* pulley1 = new WrapCylinder();
    pulley1->setName("pulley1");
    pulley1->set_radius(r);
    pulley1->set_length(0.05);

    // Add the wrap object to the body, which takes ownership of it
    ground.addWrapObject(pulley1);

    // One spring has wrap cylinder with respect to ground origin
    PathSpring* spring1 =
        new PathSpring("spring1", 1.0, 0.1, 0.01);
    spring1->updGeometryPath().
        appendNewPathPoint("origin", ground, Vec3(-off, 0, 0));
    spring1->updGeometryPath().
        appendNewPathPoint("insert", *body, Vec3(0));
    spring1->updGeometryPath().addPathWrap(*pulley1);

    model.addComponent(spring1);

    WrapCylinder* pulley2 = new WrapCylinder();
    pulley2->setName("pulley2");
    pulley2->set_radius(r);
    pulley2->set_length(0.05);

    // Add the wrap object to the body, which takes ownership of it
    bodyOffset->addWrapObject(pulley2);

    // Second spring has wrap cylinder with respect to bodyOffse origin
    PathSpring* spring2 =
        new PathSpring("spring2", 1.0, 0.1, 0.01);
    spring2->updGeometryPath().
        appendNewPathPoint("origin", ground, Vec3(-off, 0, 0));
    spring2->updGeometryPath().
        appendNewPathPoint("insert", *body, Vec3(0));
    spring2->updGeometryPath().addPathWrap(*pulley2);
    spring2->updGeometryPath().setDefaultColor(Vec3(0, 0.8, 0));

    model.addComponent(spring2);

    //model.setUseVisualizer(true);

    SimTK::State& s = model.initSystem();
    auto& coord = joint->updCoordinate();

    int nsteps = 10;
    for (int i = 0; i < nsteps; ++i) {
        
        coord.setValue(s, i*SimTK::Pi/(2*nsteps));
        model.realizeVelocity(s);

        //model.getVisualizer().show(s);

        double ma1 = spring1->computeMomentArm(s, coord);
        double ma2 = spring2->computeMomentArm(s, coord);

        ASSERT_EQUAL<double>(r, ma1, SimTK::Eps);
        ASSERT_EQUAL<double>(r, ma2, SimTK::Eps);

        double len1 = spring1->getLength(s);
        double len2 = spring2->getLength(s);

        ASSERT_EQUAL<double>(len1, len2, SimTK::Eps);
    }
}

// This additional function tests the analytical solution of the muscle length. It simulates a
// model of two vertically moving masses with in the middle a wrapping surface cylinder. A muscle
// goes from mass 1 to wrapping surface to mass 2 and is activated for a duration of 5 seconds.
// If the cylinder is not rotated (the wrapping goes straight over the x-y plane), an analytical
// solution can be computed which can be compared to the numerical one (tendonLength + fiberLength)
// TODO: extend with wrapping surfaces that go wrong now -> ellipsoids, spheres, toroid
void testWrapCylinderAnalytical(){
    // Create a new OpenSim model
    testCase tc;
    auto model = Model();
    model.setName("WrappingCylinder_AnalyticalTests");
    model.setGravity(Vec3(0, -9.80665, 0));

    // Create three bodies, two on the left and right of the cylinder, one at the floor for attachments of springs
    double bodyMass = 30.0;
    double bodySideLength = tc.BODY_SIZE;
    auto bodyInertia = bodyMass * Inertia::brick(Vec3(bodySideLength / 2.));
    auto bodyLeft = new OpenSim::Body("bodyLeft", bodyMass, Vec3(0), bodyInertia);
    auto bodyRight = new OpenSim::Body("bodyRight", bodyMass, Vec3(0), bodyInertia);
    auto bodyGround = new OpenSim::Body("bodyGround", bodyMass, Vec3(0), bodyInertia);

    model.addBody(bodyLeft);
    model.addBody(bodyRight);
    model.addBody(bodyGround);

    // Attach the pelvis to ground with a vertical slider joint, and attach the
    // pelvis, thigh, and shank bodies to each other with pin joints.
    Vec3 sliderOrientation(0, 0, SimTK::Pi / 2.);
    Vec3 bodyOffset(tc.BODY_OFFSET, 0, 0);
    auto sliderLeft = new SliderJoint("sliderLeft", model.getGround(), bodyOffset,
                                      sliderOrientation, *bodyLeft, Vec3(0), sliderOrientation);
    auto sliderRight = new SliderJoint("sliderRight", model.getGround(), -bodyOffset,
                                       sliderOrientation, *bodyRight, Vec3(0), sliderOrientation);
    auto weldGround = new WeldJoint("weldGround", model.getGround(), *bodyGround);

    // Add the joints to the model.
    model.addJoint(sliderLeft);
    model.addJoint(sliderRight);
    model.addJoint(weldGround);

    // Set the coordinate names and default values. Note that we need "auto&"
    // here so that we get a reference to the Coordinate rather than a copy.
    auto &sliderCoordLeft =
            sliderLeft->updCoordinate(SliderJoint::Coord::TranslationX);
    sliderCoordLeft.setName("yCoordSliderLeft");
    sliderCoordLeft.setDefaultValue(0.5);
    auto &sliderCoordRight =
            sliderRight->updCoordinate(SliderJoint::Coord::TranslationX);
    sliderCoordRight.setName("yCoordSliderRight");
    sliderCoordRight.setDefaultValue(0.5);

    // Create a single muscle
    double mclFmax = 4000., mclOptFibLen = 0.55, mclTendonSlackLen = 0.5,
            mclPennAng = 0.;
    auto muscle = new Thelen2003Muscle("muscle", mclFmax, mclOptFibLen,
                                       mclTendonSlackLen, mclPennAng);
    muscle->addNewPathPoint("origin", *bodyLeft, Vec3(0, bodySideLength / 2, 0));
    muscle->addNewPathPoint("insertion", *bodyRight, Vec3(0, bodySideLength / 2, 0));
    model.addForce(muscle);

    // Create the wrapping Cylinder
    auto wrappingFrame = new PhysicalOffsetFrame("wrappingFrame", model.getGround(),
                                                 SimTK::Transform(Vec3(0, tc.CYLINDER_HEIGHT, 0)));
    // Add the wrapping surface
    auto wrapSurface = new WrapCylinder();
    wrapSurface->setAllPropertiesUseDefault(true);
    wrapSurface->set_radius(tc.CYLINDER_RADIUS);
    wrapSurface->set_length(1);
    wrapSurface->set_xyz_body_rotation(Vec3(0));
    wrapSurface->set_quadrant("+y");

    wrapSurface->setName("wrapSurface");
    wrappingFrame->addWrapObject(wrapSurface);
    bodyGround->addComponent(wrappingFrame);

    // Configure the vastus muscle to wrap over the patella.
    muscle->updGeometryPath().addPathWrap(*wrapSurface);

    // Create a controller
    auto brain = new PrescribedController();
    brain->setActuators(model.updActuators());
    double t[5] = {0.0, 1.0, 2.0, 3.0, 4.0}, x[5] = {0.0, 0.3, 0.0, 0.5, 0.0};
    auto controlFunction = new PiecewiseConstantFunction(5, t, x);
    brain->prescribeControlForActuator("muscle", controlFunction);
    model.addController(brain);

    // Add table for result processing. This table will be used to compute instances of analytical solutions as well
    static const std::string sliderLPath{"/jointset/sliderLeft/yCoordSliderLeft"};
    static const std::string sliderRPath{"/jointset/sliderRight/yCoordSliderRight"};
    auto table = new TableReporter();
    table->setName("wrapping_results_table");
    table->set_report_time_interval(tc.REPORTING_INTERVAL);
    table->addToReport(model.getComponent(sliderLPath).getOutput("value"), "height slider L");
    table->addToReport(model.getComponent(sliderRPath).getOutput("value"), "height slider R");
    table->addToReport(model.getComponent("/forceset/muscle").getOutput("fiber_length"));
    table->addToReport(model.getComponent("/forceset/muscle").getOutput("tendon_length"));
    model.addComponent(table);

    SimTK::State& x0 = model.initSystem();
    // time the simulation
    clock_t ticks = clock();
    simulate(model, x0, 0.0, tc.FINAL_TIME);
    ticks = clock() - ticks;
    double runTime = (float)ticks/CLOCKS_PER_SEC;
    // Display the execution time for a simulation of 5 seconds
    cout << "Execution time: "
        << ticks << " clicks, "
        << runTime << " seconds (sim time = " << tc.FINAL_TIME << " seconds)" << endl;

    // unpack the table to analyze the results
    const auto headings = table->getTable().getColumnLabels();
    const auto leftHeight = table->getTable().getDependentColumnAtIndex(0);
    const auto rightHeight = table->getTable().getDependentColumnAtIndex(1);
    const auto fiberLength = table->getTable().getDependentColumnAtIndex(2);
    const auto tendonLength = table->getTable().getDependentColumnAtIndex(3);

    for (int i=0; i<tc.FINAL_TIME/tc.REPORTING_INTERVAL; i++){
        double lNumerical  = fiberLength[i] + tendonLength[i];

        double margin = 0.001;   // 0.1% margin allowed
        SimTK_TEST(lNumerical > (1+margin)*analyticalSolution(leftHeight[i], tc)
                   || lNumerical < (1-margin)*analyticalSolution(leftHeight[i], tc));
    }
}

double analyticalSolution(double leftHeight, const testCase& tc){
    // create shorter notation from the global parameters
    double s = tc.BODY_SIZE;
    double r = tc.CYLINDER_RADIUS;
    double x = tc.BODY_OFFSET;
    double h = tc.CYLINDER_HEIGHT;

    double lCylinder;

    if (leftHeight+s/2 >= h+r) {
        // total analytical length
        return 2*x;
    }
    else {
        double hD = h - (leftHeight + s / 2);
        // distance center cylinder and connection point
        double d = sqrt(pow(x, 2) + pow(hD, 2));
        // length from muscle connection point to tangent circle
        double lTangent = sqrt(pow(d, 2) - pow(r, 2));
        // angle connection point to horizontal
        double beta = atan(r / lTangent) + atan(hD / x);
        // height of extension of connection point to tangent circle (h + small part)
        double H = tan(beta) * x;
        // center of cylinder, horizontal line until it touches muscle
        double y = ((H - hD) / H) * x;
        // angle horizontal and perpendicular to tangent to circle
        double theta = acos(r / y);
        // length over cylinder part
        lCylinder = (SimTK::Pi - 2 * theta) * r;
        // total analytical length
        return 2 * lTangent + lCylinder;
    }
}

void simulateModelWithMusclesNoViz(const string &modelFile, double finalTime, double activation)
{
    // Create a new OpenSim model
    Model osimModel(modelFile);

    double initialTime = 0;

    // Create a prescribed controller that simply applies a function of the force
    PrescribedController actuatorController;
    actuatorController.setActuators(osimModel.updActuators());
    for (int i=0; i<actuatorController.getActuatorSet().getSize(); i++) {
        actuatorController.prescribeControlForActuator(i, new Constant(activation));
    }

    // add the controller to the model
    osimModel.addController(&actuatorController);
    osimModel.disownAllComponents(); // because PrescribedController is on stack

    osimModel.finalizeFromProperties();
    osimModel.printBasicInfo();

    // Initialize the system and get the state representing the state system
    SimTK::State& si = osimModel.initSystem();

    const Set<Muscle>& muscles = osimModel.getMuscles();
    for (int i=0; i<muscles.getSize(); i++){
        muscles[i].setActivation(si, activation); //setDisabled(si, true);
    }
    osimModel.equilibrateMuscles(si);

    simulate(osimModel, si, initialTime, finalTime);

}// end of simulateModelWithMusclesNoViz()



void simulateModelWithPassiveMuscles(const string &modelFile, double finalTime)
{
    // Create a new OpenSim model
    Model osimModel(modelFile);
    double initialTime = 0;

    // Show model visualizer
    osimModel.setUseVisualizer(true);

    // Initialize the system and get the state representing the state system
    SimTK::State& si = osimModel.initSystem();

    const Set<Muscle>& muscles = osimModel.getMuscles();
    for (int i=0; i<muscles.getSize(); ++i){
        muscles[i].setActivation(si, 0); //setDisabled(si, true);
    }
    osimModel.equilibrateMuscles(si); 

    const ModelVisualizer& modelViz = osimModel.getVisualizer();
    const Visualizer& viz = modelViz.getSimbodyVisualizer();
    viz.report(si);

    simulate(osimModel, si, initialTime, finalTime);


}// end of simulateModelWithMusclesViz()

void simulateModelWithoutMuscles(const string &modelFile, double finalTime)
{
    // Create a new OpenSim model
    Model osimModel(modelFile);
    double initialTime = 0;

    // Show model visualizer
    osimModel.setUseVisualizer(true);

    // remove all forces
    auto& forces = osimModel.updForceSet();
    for(int i = 0; i<forces.getSize(); ++i){
        forces.remove(i);
    }

    SimTK::State& si = osimModel.initSystem();

    simulate(osimModel, si, initialTime, finalTime);

}// end of simulateModelWithoutMuscles()



void simulateModelWithLigaments(const string &modelFile, double finalTime)
{
    // Create a new OpenSim model
    Model osimModel(modelFile);
    double initialTime = 0;

    // Show model visualizer
    osimModel.setUseVisualizer(true);

    // TODO replace muscles with ligaments
//    Set<OpenSim::Force> &forces = osimModel.updForceSet();
//    for(int i = 0; i<forces.getSize(); ++i){
//        forces.remove(&forces[i]);
//    }

    SimTK::State& si = osimModel.initSystem();

    simulate(osimModel, si, initialTime, finalTime);

}// end of simulateModelWithoutMuscles()




void simulateModelWithCables(const string &modelFile, double finalTime)
{
    // Create a new OpenSim model
    Model osimModel(modelFile);
    double initialTime = 0;

    SimTK::State& si = osimModel.initSystem();

    Array<GeometryPath*> paths;
    Array<String> pathNames;
    int numLigaments = 0, numMuscles = 0;
    for (int i = 0; i < osimModel.getForceSet().getSize(); ++i) {
        Ligament* lig = dynamic_cast<Ligament*>(&osimModel.getForceSet()[i]);
        if (lig != 0) {
            numLigaments++;
            paths.append(&lig->updGeometryPath());
            pathNames.append(lig->getName());
            continue;
        }

        Muscle* mus = dynamic_cast<Muscle*>(&osimModel.getForceSet()[i]);
        if (mus != 0) {
            numMuscles++;
            paths.append(&mus->updGeometryPath());
            pathNames.append(mus->getName());
            cout << mus->getName() << ": " << mus->getGeometryPath().getWrapSet().getSize() << endl;
            continue;
        }
    }

    cout << "num ligaments = " << numLigaments << endl;
    cout << "num muscles = " << numMuscles << endl;

    Array<CableInfo> cableInfos;
    std::set<String> wrapObjectsInUse;
    for (int i = 0; i < paths.getSize(); ++i) {
        CableInfo cableInfo;

        GeometryPath* geomPath = paths[i];
        const PathWrapSet& wrapSet = geomPath->getWrapSet();
        const PathPointSet& viaSet = geomPath->getPathPointSet();
        Array<AbstractPathPoint*> activePathPoints = geomPath->getCurrentPath(si);

        AbstractPathPoint* orgPoint = &viaSet[0];
        AbstractPathPoint* insPoint = &viaSet[viaSet.getSize()-1];

        cableInfo.orgBodyName = orgPoint->getParentFrame().getName();
        cableInfo.orgLoc = orgPoint->getLocation(si);
        cableInfo.insBodyName =  insPoint->getParentFrame().getName();
        cableInfo.insLoc = insPoint->getLocation(si);

        int numVias = 0, numSurfs = 0;
        for (int j = 0; j < activePathPoints.getSize(); ++j) {
            AbstractPathPoint* pp = activePathPoints[j];
            cout << "pp" << j << " = " << pp->getName() << " @ "
                 << pp->getParentFrame().getName() << ", loc = "
                 << pp->getLocation(si) << endl;
        }

        for (int j = 0; j < viaSet.getSize(); ++j) {
            AbstractPathPoint* pp = &viaSet[j];
            cout << "mp" << j << " = " << pp->getName() << " @ "
                 << pp->getParentFrame().getName() << ", loc = "
                 << pp->getLocation(si) << endl;
        }

        // add vias to cableInfo
        for (int j = 1; j < viaSet.getSize()-1; ++j) { // skip first (origin) and last (insertion) points
            const AbstractPathPoint& pp = viaSet[j];
            ObstacleInfo obs;
            obs.bodyName = pp.getParentFrame().getName();
            obs.X_BS.setP(pp.getLocation(si));
            obs.wrapObjectPtr=NULL;
            obs.P_S.setToNaN(); // not used for viapoint
            obs.Q_S.setToNaN(); // not used for viapoint
            obs.isVia = true;
            obs.isActive = true;
            cableInfo.obstacles.append(obs);
            numVias++;
        }

        Array<ObstacleInfo*> wrapObs;
        // add surfaces to cableInfo
        for (int j = 0; j < wrapSet.getSize(); j++) {
            const PathWrap& wrap = wrapSet[j];
            const WrapObject* wrapObj = wrap.getWrapObject();
            ObstacleInfo obs;
            obs.wrapObjectPtr = wrapObj;
            obs.bodyName = wrapObj->getFrame().getName();
            obs.X_BS = wrapObj->getTransform();
            obs.isVia = false;
            // initially assume inactive
            obs.isActive = false;
            obs.P_S.setToNaN();
            obs.Q_S.setToNaN();
            wrapObs.append(&obs);

            wrapObjectsInUse.insert(wrapObj->getName());
            numSurfs++;
        }

        // find active wrap surfaces and insert to obstacles list in order
        int obsIdx = 0;
        int viaPtIdx = 1; // skip first (origin) point
        for (int j = 1; j < activePathPoints.getSize()-1; ++j) { // skip first (origin) and last (insertion)
            AbstractPathPoint* pp = activePathPoints[j];
            if (*pp == viaSet[viaPtIdx]) {
                viaPtIdx++;
                obsIdx++;
                continue;
            }
            else { // next two path points should be a wrap point
                for (int k = 0; k < wrapSet.getSize(); ++k) {
                    const Vec3& wrapStartPointLoc = wrapSet[k].getPreviousWrap().r1;
                    if (!wrapStartPointLoc.isInf() && pp->getLocation(si).isNumericallyEqual(wrapStartPointLoc)) {
                        ObstacleInfo* obs = wrapObs[k];
                        obs->isActive = true;
                        // pp and next pp are wrap points
                        Transform X_SB = obs->X_BS.invert();
                        AbstractPathPoint* pp_next = activePathPoints[++i]; // increment to next pp
                        obs->P_S = X_SB*pp->getLocation(si);
                        obs->Q_S = X_SB*pp_next->getLocation(si);
                        cableInfo.obstacles.insert(obsIdx++, *obs);
                        break;
                    }
                }
            }
        }

        // append inactive wrap surfaces to obstacle list
        // TODO find order of inactive wrap surfaces relative to via points
        for (int j = 0; j < wrapObs.getSize(); ++j) {
            ObstacleInfo* obs = wrapObs[j];
            if (!obs->isActive) {
                cableInfo.obstacles.append(*obs);
            }
        }

        cableInfos.append(cableInfo);

        cout << pathNames[i] << " path, num pts " << activePathPoints.getSize()
            << ", num vias = " << numVias << ", num surfs = " << numSurfs << endl;


        // debugging

        cout << "org" << " = " << orgPoint->getName() << " @ " << 
            orgPoint->getParentFrame().getName() << ", loc = " 
            << orgPoint->getLocation(si) << endl;
        cout << "ins" << " = " << insPoint->getName() << " @ " <<
            insPoint->getParentFrame().getName() << ", loc = " 
            << insPoint->getLocation(si) << endl;

        for (int j = 0; j < cableInfo.obstacles.getSize(); ++j) {
            ObstacleInfo oi = cableInfo.obstacles[j];
            cout << "ObsInfo " << j << ": via=" << oi.isVia << ", bodyName=" << oi.bodyName << ", loc=" << oi.X_BS.p() << ", P=" << oi.P_S << ", Q=" << oi.Q_S << endl;
        }


//        for (int j = 0; j < wrapSet.getSize(); ++j) {
//            cout << "wrap object " << j << " name = " << wrapSet[j].getName() << endl;
//            cout << "wrap point 0 = " << wrapSet[j].getWrapPoint(0).getLocation() << endl;
//            cout << "wrap point 1 = " << wrapSet[j].getWrapPoint(1).getLocation() << endl;
//            const WrapResult& wr = wrapSet[j].getPreviousWrap();
//            cout << "wrap result r1 = " << wr.r1 << endl;
//            cout << "wrap result r2 = " << wr.r2 << endl;
//            cout << "wrap result startpt = " << wr.startPoint << endl;
//            cout << "wrap result endpt = " << wr.endPoint << endl;
//
////            for (int j = 0; j < wr.wrap_pts.getSize(); j++) {
////                cout << "wrap result pt[" << j << "] = " << wr.wrap_pts[j] << endl;
////            }
//        }
    }


    cout << "num forces before removal = " << osimModel.getForceSet().getSize() << endl;

    // remove all OpenSim forces
    osimModel.updForceSet().clearAndDestroy();

    // remove all unused wrap objects
    cout << "wraps in use:" << endl;
    set<String>::iterator it;
    for ( it=wrapObjectsInUse.begin() ; it != wrapObjectsInUse.end(); it++ ) {
        cout << "    " << *it << endl;
    }

    for (int i = 0; i < osimModel.getBodySet().getSize(); ++i) {
        Array<WrapObject*> toRemove;
        const WrapObjectSet& wraps = osimModel.getBodySet()[i].getWrapObjectSet();
        // XXX hack to remove wrap objects, which is not supported in the API
        WrapObjectSet *mutableWraps = const_cast<WrapObjectSet *>(&wraps);
        cout << "culling wrap objs from " << osimModel.getBodySet()[i].getName() << ", num wraps = " << wraps.getSize() << endl;
//        mutableWraps->clearAndDestroy();
        mutableWraps->setMemoryOwner(true);
        for (int j = 0; j < wraps.getSize(); ++j) {
            if (wrapObjectsInUse.find(wraps[j].getName()) == wrapObjectsInUse.end()) { // does not contain
                toRemove.append(&wraps[j]);
            }
        }

        cout << "wraps to remove:" << endl;
        for (int j = 0; j < toRemove.getSize(); ++j) {
            cout << "    " << toRemove[j]->getName() << endl;
        }

//        cout << "indices to remove = " << toRemove << endl;
        for (int j = 0; j < toRemove.getSize(); ++j) {
            mutableWraps->remove(toRemove[j]);
        }
        cout << "finished cull from " << osimModel.getBodySet()[i].getName() << ", num wraps = " << wraps.getSize() << endl;

    }

//    cout << "counting wraps in testWrapping..." << endl;
//    for (int i = 0; i < osimModel.getBodySet().getSize(); i++) {
//           const OpenSim::Body& body = osimModel.getBodySet()[i];
//           const WrapObjectSet& wrapObjects = body.getWrapObjectSet();
//           cout << body.getName() << ": wrap count = " << wrapObjects.getSize() << endl;
//    }

    // Show model visualizer
    osimModel.setUseVisualizer(true);

    // Initialize the system and get the state representing the state system
//    SimTK::State& s = osimModel.initSystem();
    osimModel.buildSystem();

    const ModelVisualizer& modelViz = osimModel.getVisualizer();
    const Visualizer& viz = modelViz.getSimbodyVisualizer();
    viz.setBackgroundColor(Vec3(1,1,1));

    cout << "num forces after removal = " << osimModel.getForceSet().getSize() << endl;


    MultibodySystem& system = osimModel.updMultibodySystem();
    GeneralForceSubsystem& forceSubsystem = osimModel.updForceSubsystem();
    // const SimbodyMatterSubsystem& matter = osimModel.getMatterSubsystem();
    // const SimbodyEngine& engine = osimModel.getSimbodyEngine();

    // add cable system
    CableTrackerSubsystem cables(system);

    for (int i = 0; i < cableInfos.getSize(); ++i) {
        const CableInfo& cableInfo = cableInfos[i];

        const OpenSim::Body& orgBody = osimModel.getBodySet().get(cableInfo.orgBodyName);
        const OpenSim::Body& insBody = osimModel.getBodySet().get(cableInfo.insBodyName);
        const MobilizedBody& orgMobBody = orgBody.getMobilizedBody();
        const MobilizedBody& insMobBody = insBody.getMobilizedBody();
//        cout << "origin mob idx = " << orgBody.getIndex() << endl;
//        cout << "insertion mob idx = " << insBody.getIndex() << endl;


        CablePath path(cables, orgMobBody, cableInfo.orgLoc,   // origin
                            insMobBody, cableInfo.insLoc);  // termination

        // Add obstacles
        for (int j = 0; j < cableInfo.obstacles.getSize(); ++j) {
            ObstacleInfo oi = cableInfo.obstacles[j];
            const OpenSim::Body& osBody = osimModel.getBodySet().get(oi.bodyName);
            const MobilizedBody& mobBody = osBody.getMobilizedBody();

            if (oi.isVia) {
                CableObstacle::ViaPoint via(path, mobBody, oi.X_BS.p());
            }
            else {
                const WrapSphere* wrapSphere = dynamic_cast<const WrapSphere*>(oi.wrapObjectPtr);
                if (wrapSphere != 0) {
                    CableObstacle::Surface surf(path, mobBody,
                        oi.X_BS, SimTK::ContactGeometry::Sphere(wrapSphere->getRadius())); // along y
                    if (oi.isActive)
                        surf.setContactPointHints(oi.P_S, oi.Q_S);
                    else
                        surf.setDisabledByDefault(true);
                    continue;
                }

                const WrapCylinder* wrapCyl = dynamic_cast<const WrapCylinder*>(oi.wrapObjectPtr);
                if (wrapCyl != 0) {
                    CableObstacle::Surface surf(path, mobBody,
                        oi.X_BS, SimTK::ContactGeometry::Cylinder(wrapCyl->get_radius())); // along y
                    if (oi.isActive)
                        surf.setContactPointHints(oi.P_S, oi.Q_S);
                    else
                        surf.setDisabledByDefault(true);
                    continue;
                }

                const WrapEllipsoid* wrapEllip = dynamic_cast<const WrapEllipsoid*>(oi.wrapObjectPtr);
                if (wrapEllip != 0) {
                    CableObstacle::Surface surf(path, mobBody,
                        oi.X_BS, SimTK::ContactGeometry::Ellipsoid(wrapEllip->getRadii())); // along y
                    if (oi.isActive)
                        surf.setContactPointHints(oi.P_S, oi.Q_S);
                    else
                        surf.setDisabledByDefault(true);
                    continue;
                }

                const WrapTorus* wrapTorus = dynamic_cast<const WrapTorus*>(oi.wrapObjectPtr);
                if (wrapTorus != 0) {
                    Real tubeRadius = (wrapTorus->getOuterRadius()-wrapTorus->getInnerRadius())/2.0;
                    Real torusRadius = wrapTorus->getOuterRadius()-tubeRadius;
                    CableObstacle::Surface surf(path, mobBody,
                        oi.X_BS, SimTK::ContactGeometry::Torus(torusRadius, tubeRadius)); // along y
                    if (oi.isActive)
                        surf.setContactPointHints(oi.P_S, oi.Q_S);
                    else
                        surf.setDisabledByDefault(true);
                    continue;
                }

                cout << "unknown wrap object [" << oi.wrapObjectPtr->getName() << "], adding via point at origin" << endl;
                CableObstacle::ViaPoint via(path, mobBody, oi.X_BS.p());

            }
        }

        MyCableSpring cable(forceSubsystem, path, 100., 3.5, 0*0.1);
        //    system.addEventReporter(new ShowStuff(system, cable, 0.1*0.1));

    }

    system.addEventReporter(new ShowStuff(system, 0.1*0.1));

    SimTK::State& state = osimModel.initializeState();
    viz.report(state);
    for (int i = 0; i < cables.getNumCablePaths(); ++i) {
        const CablePath& path = cables.getCablePath((CablePathIndex)i);
        cout << "path " << i << " length = " << path.getCableLength(state) << endl;
    }

//    cout << "Hit ENTER ...";
//    cout.flush();
//    getchar();

    // Simulate it.
    saveStates.clear(); saveStates.reserve(2000);

    ShowStuff::showHeading(cout);

    simulate(osimModel, state, initialTime, finalTime);

}// end of simulateModelWithCables()

void simulate(Model& osimModel, State& si, double initialTime, double finalTime)
{
    // Dump model back out; no automated test provided here though.
    // osimModel.print(osimModel.getName() + "_out.osim");

    // Create the Manager for the simulation.
    const double accuracy = 1.0e-4;
    Manager manager(osimModel);
    manager.setIntegratorAccuracy(accuracy);

    // Integrate from initial time to final time
    si.setTime(initialTime);
    cout << "\nIntegrating from " << initialTime << " to " << finalTime << endl;

    const double start = SimTK::realTime();
    manager.initialize(si);
    manager.integrate(finalTime);
    cout << "simulation time = " << SimTK::realTime()-start
         << " seconds (wallclock time)\n" << endl;

    auto& integrator = manager.getIntegrator();
    cout << "integrator iterations = " << integrator.getNumStepsTaken() << endl;

    // Save the simulation results
    Storage states(manager.getStateStorage());
    states.print(osimModel.getName()+"_states.sto");
    osimModel.updSimbodyEngine().convertRadiansToDegrees(states);
    states.setWriteSIMMHeader(true);
    states.print(osimModel.getName()+"_states_degrees.mot");
} // end of simulate()

// In XMLDocument version 30515, we converted VisibleObject, color and
// display_preference properties to Appearance properties.
void testWrapObjectUpdateFromXMLNode30515() {
    XMLDocument doc("testWrapObject_updateFromXMLNode30515.osim");

    // Make sure this test is not weakened by the model in the repository being
    // updated.
    SimTK_TEST(doc.getDocumentVersion() == 20302);
    Model model("testWrapObject_updateFromXMLNode30515.osim");
    model.print("testWrapObject_updateFromXMLNode30515_updated.osim");
    const auto& wrapObjSet = model.getGround().getWrapObjectSet();

    // WrapSphere has:
    //   display_preference = 1
    //   color = default
    //   VisibleObject display_preference = 0
    {
        const auto& sphere = wrapObjSet.get("wrapsphere");
        SimTK_TEST(!sphere.get_Appearance().get_visible());
        SimTK_TEST_EQ(sphere.get_Appearance().get_color(), SimTK::Cyan);
        SimTK_TEST(sphere.get_Appearance().get_representation() == 
                VisualRepresentation::DrawPoints /* == 1 */);
    }

    // WrapCylinder has:
    //   display_preference = 0
    //   color = default
    //   VisibleObject display_preference = 4 
    {
        const auto& cyl = wrapObjSet.get("wrapcylinder");
        // The outer display_preference overrides the inner one.
        SimTK_TEST(!cyl.get_Appearance().get_visible());
        SimTK_TEST_EQ(cyl.get_Appearance().get_color(), SimTK::Cyan);
        SimTK_TEST(cyl.get_Appearance().get_representation() == 
                VisualRepresentation::DrawSurface /* == 3 */);
    }

    // WrapEllipsoid has:
    //   display_preference = 2
    //   color = 1 0.5 0
    //   VisibleObject display_preference = 3
    {
        const auto& ellipsoid = wrapObjSet.get("wrapellipsoid");
        SimTK_TEST(ellipsoid.get_Appearance().get_visible());
        SimTK_TEST_EQ(ellipsoid.get_Appearance().get_color(), Vec3(1, 0.5, 0));
        // Outer display_preference overrides inner one.
        SimTK_TEST(ellipsoid.get_Appearance().get_representation() == 
                VisualRepresentation::DrawWireframe /* == 2 */);
    }
}





