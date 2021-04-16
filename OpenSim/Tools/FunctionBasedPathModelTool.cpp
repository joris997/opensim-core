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
    // skip exe name
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
                case 0:  // MODEL
                    sourceModelPath = arg;
                    break;
                case 1:  // OUTPUT_NAME
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

            // handle flagged args (e.g. --arg)
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

        assert(argc == 0);
        assert(sourceModelPath != nullptr);
        assert(outputName != nullptr);
    }

    Model sourceModel{sourceModelPath};
    Model outputModel{sourceModel};

    sourceModel.finalizeConnections();
    sourceModel.finalizeFromProperties();
    sourceModel.initSystem();
    if (verbose) {
        sourceModel.printSubcomponentInfo();
    }

    // struct that holds how a PBP in the source maps onto an actuator in the
    // destination
    struct PBPtoActuatorMapping final {
        PointBasedPath& sourcePBP;
        PathActuator& outputActuator;

        PBPtoActuatorMapping(PointBasedPath& sourcePBP_,
                             PathActuator& outputActuator_) :
            sourcePBP{sourcePBP_},
            outputActuator{outputActuator_} {
        }
    };

    // find PBPs in the source and figure out how they map onto `PathActuator`s
    // in the destination
    //
    // (this is because `PathActuator`s are the "owners" of `GeometryPath`s in
    //  most models)
    std::vector<PBPtoActuatorMapping> mappings;
    for (auto& pa : sourceModel.updComponentList<PathActuator>()) {

        PointBasedPath* pbp = dynamic_cast<PointBasedPath*>(&pa.updGeometryPath());

        // if the actuator doesn't use a PBP, ignore it
        if (!pbp) {
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

        mappings.emplace_back(*pbp, *paDest);
    }

    // for each `PathActuator` that uses a PBP, create an equivalent
    // `FunctionBasedPath` (FBP) by fitting a function against the PBP and
    // replace the PBP owned by the destination's `PathActuator` with the FBP
    int i = 1;
    for (const auto& mapping : mappings) {
        // create an FBP in-memory
        FunctionBasedPath fbp = FunctionBasedPath::fromPointBasedPath(sourceModel, mapping.sourcePBP);

        // write the FBP's data to the filesystem
        std::string dataFilename = [&outputName, &i]() {
            std::stringstream ss;
            ss << outputName << "_DATA_" << i <<  ".txt";
            return std::move(ss).str();
        }();

        std::ofstream dataFile{dataFilename};

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

        // update the FBP to refer to the data file
        fbp.setDataPath(dataFilename);

        // assign the FBP over the destination's PBP
        mapping.outputActuator.updProperty_GeometryPath().setValue(fbp);

        ++i;
    }

    // the output model is now the same as the source model, but each PBP in
    // its `PathActuator`s has been replaced with an FBP. Perform any final
    // model-level fixups and save the output model.

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
