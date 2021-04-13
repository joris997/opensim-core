#include <vector>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstddef>
#include <regex>
#include <OpenSim/OpenSim.h>

using namespace OpenSim;

static const char HELP[] =
        "This tool converts a model with classic GeometryPaths or \
        PointBasedPaths (the path of the PathActuator is described \
        by discrete sets of points) to a model with FunctionBasedPaths \
        \
        To that end, when a model is provided, this tool creates a new \
        .osim file according to the convention -originalName-FBP.osim \
        and a .txt file containing the sampled interpolation grid data \
        \
        usage: \
            FunctionBasedPathModelTool [--help] pathToModel/model.osim \
        e.g.: \
            FunctionBasedPathModelTool /home/user/models/Arm26.osim";

int main(int argc, char **argv){
    // skip tool name
    --argc;
    ++argv;

    // parse input arguments
    std::string pbpModelString;
    if (argc == 0){
        std::cout << "Provide an argument ('--help')" << std::endl;
        return 0;
    } else if (argc > 1) {
        std::cout << "Provided too many arguments ('--help')" << std::endl;
        return 0;
    } else {
        char const* arg = argv[0];
        if (!strcmp(argument,"--help")){
            std::cout << HELP;
            return 0;
        } else {
            pbpModelString = argv[0];
        }
    }
    Model pbpModel(pbpModelString);
    Model* fbpModel = pbpModel.clone();

    pbpModel.finalizeConnections();
    pbpModel.finalizeFromProperties();
    pbpModel.initSystem();
    pbpModel.printSubcomponentInfo();

    // Create folder for txt file data and new osim file
    std::string modelPath = pbpModelString;
    std::size_t found = modelPath.find_last_of("/\\");
    std::string modelName = modelPath.substr(found+1);
    modelName = modelName.substr(0,modelName.find(".",0));

    if (mkdir(modelName.c_str(),0777) == -1){
        std::cout << "Folder for model; " << modelName << " already created";
        std::cout << "so I'll just use that one and overwrite everything"
                  << std::endl;
    } else {
        std::cout << "Folder for model; " << modelName << " created"
                  << std::endl;
    }

    // Converting
    std::vector<FunctionBasedPath> fbps;
    std::vector<PointBasedPath> pbps;
    std::ofstream printFile;
    int id = 0;

    for(PathActuator& pa :
        pbpModel.updComponentList<PathActuator>()){
        // const reference to a geometrypath
        auto &pbp = pa.getGeometryPath();

        // dynamic cast to check if its a pointbasedpath
        const PointBasedPath* pbpp = dynamic_cast<const PointBasedPath*>(&pbp);
        pbps.push_back(*pbpp);
        if(pbpp != nullptr){
            std::cout << "pa: " << pa.getName() << std::endl;

            FunctionBasedPath fbp(pbpModel,*pbpp,id);

            printFile.open(modelName+"/FBP"+std::to_string(id)+".txt");
            fbp.printContent(printFile);

            fbps.push_back(fbp);
            id++;
        }
    }

    int cnt = 0;
    const PointBasedPath* pbpp;
    std::vector<PathActuator*> pap;
    for(PathActuator& pa :
        fbpModel->updComponentList<PathActuator>()){
        pap.push_back(&pa);
    }
    for (unsigned i=0; i<pap.size(); i++){
        auto &pbp = pap[i]->getGeometryPath();
        pbpp = dynamic_cast<const PointBasedPath*>(&pbp);
        if(pbpp != nullptr){
            pap[i]->updProperty_GeometryPath().setValue(fbps[cnt]);
            cnt++;
        }
    }
    fbpModel->finalizeConnections();
    fbpModel->finalizeFromProperties();
    fbpModel->initSystem();
    fbpModel->printSubcomponentInfo();

    // Print model to file
    fbpModel->print(modelName+"FBP.osim");

    return 0;
}
