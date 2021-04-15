#include <OpenSim/OpenSim.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstddef>
#include <sstream>

using namespace OpenSim;

static const char usage[] = "FunctionBasedPathModelTool [--help] MODEL OUTPUT_NAME";
static const char help_text[] =
R"(
ARGS

    --help           print this help
    -v, --verbose    print verbose output

DESCRIPTION

    This tool tries to convert an .osim file (MODEL) that contains
    `GeometryPath`s and `PointBasedPath`s, into a functionally equivalent model
    file, written to ${OUTPUT_NAME}.osim, that contains `FunctionBasedPath`s.

    The main reason to do this is performance. Computing a path at each timestep
    in a simulation can be computationally expensive - especially if the path is
    complex (e.g. contains complex wrapping geometry). Pre-computing the
    path's topology ahead of time into a function-based lookup *may* speed up
    simulations that spend a significant amount of time computing paths.

    This tool effectively:

      - Reads MODEL (the source model)
      - Searches for `GeometryPath`s and `PointBasedPath`s in the source model
      - Tries to parameterize those paths against the source model's coordinates
        (e.g. joint coordinates), to produce an n-dimensional Bezier curve fit
        of those paths
      - Saves the fit data to ${OUTPUT_NAME}_DATA_${i}.txt, where `i` is an
        arbirary ID that links the `FunctionBasedPath` in the output .osim file
        to the Bezier fit's data
      - Updates the source model to contain `FunctionBasedPath`s
      - Writes the updated model to `${OUTPUT_NAME}.osim`

EXAMPLE USAGE

    FunctionBasedPathModelTool RajagopalModel.osim RajagopalModel_Precomputed.osim
)";

int main(int argc, char **argv) {
    // skip "FunctionBasedPathModelTool"
    --argc;
    ++argv;

    bool verbose = false;
    char const* sourceModelPath = nullptr;
    char const* outputName = nullptr;

    // parse CLI args
    {
        int nUnnamed = 0;
        while (argc) {
            const char* arg = argv[0];

            // handle unnamed args
            if (arg[0] != '-') {
                switch (nUnnamed) {
                case 0:
                    sourceModelPath = arg;
                    break;
                case 1:
                    outputName = arg;
                    break;
                default:
                    std::cerr << "FunctionBasedPathModelTool: error: too many arguments: should only supply 'MODEL' and 'OUTPUT_NAME'\n\nUSAGE: "
                              << usage
                              << std::endl;
                    return -1;
                }
                ++nUnnamed;
                --argc;
                ++argv;
                continue;
            }

            // flagged args
            if (std::strcmp(arg, "--help") == 0) {
                std::cout << "usage: " << usage << '\n' << help_text << std::endl;
                return 0;
            } else if (std::strcmp(arg, "--verbose") == 0 || std::strcmp(arg, "-v") == 0) {
                verbose = true;
                --argc;
                ++argv;
            } else {
                std::cerr << "FunctionBasedPathModelTool: error: unknown argument '"
                          << arg
                          << "': see --help for usage info";
                return -1;
            }
        }

        if (nUnnamed != 2) {
            std::cerr << "FunctionBasedPathModelTool: error: too few arguments supplied\n\n"
                      << usage
                      << std::endl;
            return -1;
        }
    }

    Model sourceModel{sourceModelPath};
    Model outputModel{sourceModel};

    sourceModel.finalizeConnections();
    sourceModel.finalizeFromProperties();
    sourceModel.initSystem();
    if (verbose) {
        sourceModel.printSubcomponentInfo();
    }

    // struct for holding source + dest actuators affected by this process
    struct PathActuatorWithPointBasedPath final {
        PathActuator& inSource;
        PathActuator& inOutput;

        PathActuatorWithPointBasedPath(PathActuator& inSource_,
                                       PathActuator& inOutput_) :
            inSource{inSource_},
            inOutput{inOutput_} {
        }
    };

    // collect actuators with `PointBasedPath`s (PBPs)
    std::vector<PathActuatorWithPointBasedPath> actuators;
    for (auto& pa : sourceModel.updComponentList<PathActuator>()) {

        // if the actuator doesn't use a PBP, ignore it
        if (!dynamic_cast<PointBasedPath*>(&pa.updGeometryPath())) {
            continue;
        }

        // otherwise, find the equivalent path in the destination
        Component* c = const_cast<Component*>(outputModel.findComponent(pa.getAbsolutePath()));

        if (!c) {
            std::stringstream err;
            err << "could not find '" << pa.getAbsolutePathString() << "' in the output model: this is a programming error";
            throw std::runtime_error{std::move(err).str()};
        }

        PathActuator* paDest = dynamic_cast<PathActuator*>(c);

        if (!paDest) {
            std::stringstream err;
            err << "the component '" << pa.getAbsolutePathString() << "' has a class of '" << pa.getConcreteClassName() << "' in the source model but a class of '" << c->getConcreteClassName() << "' in the destination model: this shouldn't happen";
            throw std::runtime_error{std::move(err).str()};
        }

        actuators.emplace_back(pa, *paDest);
    }

    // for each actuator with a PBP, create an equivalent `FunctionBasedPath`
    // and handle saving the fits etc.
    int i = 1;
    for (const auto& actuator : actuators) {
        // this should be a safe cast, because it is checked above
        const auto& pbp =
            dynamic_cast<const PointBasedPath&>(actuator.inSource.getGeometryPath());

        // create an FBP in-memory
        FunctionBasedPath fbp = FunctionBasedPath::fromPointBasedPath(sourceModel, pbp);

        // write the FBP's data to the filesystem
        std::string dataFilename;
        {
            std::stringstream ss;
            ss << outputName << "_DATA_" << i <<  ".txt";
            dataFilename = std::move(ss).str();
        }

        std::fstream dataFile{dataFilename};

        if (!dataFile) {
            std::stringstream msg;
            msg << "error: could not open a `FunctionBasedPath`'s data file at: " << dataFilename;
            throw std::runtime_error{std::move(msg).str()};
        }

        fbp.printContent(dataFile);

        if (!dataFile) {
            std::stringstream msg;
            msg << "error: error occurred after writing `FunctionBasedPath`s data to" << dataFilename;
            throw std::runtime_error{std::move(msg).str()};
        }

        fbp.setDataPath(dataFilename);

        // assign the FBP to the output model
        actuator.inOutput.updProperty_GeometryPath().setValue(fbp);

        ++i;
    }

    // FBPs substitued: perform any final steps and write the output model

    outputModel.finalizeFromProperties();
    outputModel.initSystem();
    outputModel.print(std::string{outputName} + ".osim");

    if (verbose) {
        std::cerr << "--- interpolation complete ---\n\n"
                  << "model before:\n";
        sourceModel.printSubcomponentInfo();
        std::cerr << "\nmodel after:\n";
        outputModel.printSubcomponentInfo();
    }

    return 0;
}
