#include <OpenSim/OpenSim.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstddef>
#include <sstream>

using namespace OpenSim;
using namespace std;

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

int main(int argc, char **argv)
{
    try {
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
        FunctionBasedPathConversionTool tool(sourceModelPath,outputName);
        bool success = tool.run();
        if (success){
            return 0;
        } else {
            return 1;
        }

    } catch(const std::exception& x) {
        log_error("Exception in ID: {}", x.what());
        return -1;
    }
    return 0;
}
