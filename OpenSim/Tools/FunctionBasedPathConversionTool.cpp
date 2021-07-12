#include <OpenSim/OpenSim.h>
#include "FunctionBasedPathConversionTool.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <cstddef>
#include <sstream>

using namespace OpenSim;

FunctionBasedPathConversionTool::~FunctionBasedPathConversionTool()
{
//    delete &_modelPath;
//    delete &_newModelName;
}

FunctionBasedPathConversionTool::FunctionBasedPathConversionTool()
    : _modelPath(""), _newModelName("")
{

}

FunctionBasedPathConversionTool::FunctionBasedPathConversionTool(
        const std::string modelPath, const std::string newModelName)
    : _modelPath(modelPath), _newModelName(newModelName)
{

}

bool FunctionBasedPathConversionTool::run(){
    try{
        Model sourceModel{_modelPath};
        Model outputModel{sourceModel};

        sourceModel.finalizeConnections();
        sourceModel.finalizeFromProperties();
        sourceModel.initSystem();

        outputModel.finalizeConnections();
        outputModel.finalizeFromProperties();

        bool verbose = false;
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
                throw std::runtime_error{move(err).str()};
            }

            PathActuator* paDest = dynamic_cast<PathActuator*>(c);

            if (!paDest) {
                std::stringstream err;
                err << "the component '" << pa.getAbsolutePathString() << "' has a class of '" << pa.getConcreteClassName() << "' in the source model but a class of '" << c->getConcreteClassName() << "' in the destination model: this shouldn't happen";
                throw std::runtime_error{move(err).str()};
            }

            mappings.emplace_back(*pbp, *paDest);
        }

        // for each `PathActuator` that uses a PBP, create an equivalent
        // `FunctionBasedPath` (FBP) by fitting a function against the PBP and
        // replace the PBP owned by the destination's `PathActuator` with the FBP
        int i = 1;
        std::string newModelNameLocal = _newModelName;
        for (const auto& mapping : mappings) {
            // create an FBP in-memory
            FunctionBasedPath fbp = FunctionBasedPath::fromPointBasedPath(sourceModel, mapping.sourcePBP);

            // write the FBP's data to the filesystem
            std::string dataFilename = [newModelNameLocal, &i]() {
                std::stringstream ss;
                ss << newModelNameLocal << "_DATA_" << i <<  ".txt";
                return move(ss).str();
            }();

            std::ofstream dataFile{dataFilename};

            if (!dataFile) {
                std::stringstream msg;
                msg << "error: could not open a `FunctionBasedPath`'s data file at: " << dataFilename;
                throw std::runtime_error{move(msg).str()};
            }

            fbp.printContent(dataFile);

            if (!dataFile) {
                std::stringstream msg;
                msg << "error: error occurred after writing `FunctionBasedPath`s data to" << dataFilename;
                throw std::runtime_error{move(msg).str()};
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
        outputModel.finalizeConnections();
        outputModel.initSystem();
        outputModel.print(std::string{_newModelName} + ".osim");

        if (verbose) {
            std::cerr << "--- interpolation complete ---\n\n"
                      << "model before:\n";
            sourceModel.printSubcomponentInfo();
            std::cerr << "\nmodel after:\n";
            outputModel.printSubcomponentInfo();
        }
    } catch(const std::exception& x) {
        log_error("Exception in FunctionBasedPathConversionTool: run", x.what());
        return false;
    }

    return true;
}
