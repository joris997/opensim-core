#ifndef __FunctionBasedPathConversionTool_h__
#define __FunctionBasedPathConversionTool_h__
/* -------------------------------------------------------------------------- *
 *                   OpenSim:  FunctionBasedPathConversionTool.h              *
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
#include <OpenSim/Simulation/Model/AbstractTool.h>

#include "osimToolsDLL.h"

#ifdef SWIG
    #ifdef OSIMTOOLS_API
        #undef OSIMTOOLS_API
        #define OSIMTOOLS_API
    #endif
#endif

namespace OpenSim {
//=============================================================================
//=============================================================================
/**
 * A concrete tool transforming a pointbasedpath model to a functionbasedpath
 * model
 *
 * @author Joris Verhagen
 * @version 1.0
 */
class OSIMTOOLS_API FunctionBasedPathConversionTool: public AbstractTool {
OpenSim_DECLARE_CONCRETE_OBJECT(FunctionBasedPathConversionTool, AbstractTool);

//=============================================================================
// MEMBER VARIABLES
//=============================================================================
protected:
    std::string _modelPath;
    std::string _newModelName;

//=============================================================================
// METHODS
//=============================================================================
    //--------------------------------------------------------------------------
    // CONSTRUCTION
    //--------------------------------------------------------------------------
public:
    virtual ~FunctionBasedPathConversionTool();
    FunctionBasedPathConversionTool();
    FunctionBasedPathConversionTool(const std::string modelPath, const std::string newModelName);

//    void generateModelAndPrint();

    //--------------------------------------------------------------------------
    // OPERATORS
    //--------------------------------------------------------------------------
public:
//#ifndef SWIG
//    FunctionBasedPathConversionTool&
//        operator=(const FunctionBasedPathConversionTool &aFunctionBasedPathConversionTool);
//#endif
//    void updateFromXMLNode(SimTK::Xml::Element& aNode, int versionNumber=-1) override;

    //--------------------------------------------------------------------------
    // GET AND SET
    //--------------------------------------------------------------------------
    void setModelPath(const std::string modelPath) {_modelPath = modelPath;}
    void setNewModelName(const std::string newModelName) {_newModelName = newModelName;}
    std::string getModelPath() {return _modelPath;}
    std::string getNewModelName() {return _newModelName;}

    //--------------------------------------------------------------------------
    // INTERFACE
    //--------------------------------------------------------------------------
    bool run() override SWIG_DECLARE_EXCEPTION;

//=============================================================================
};  // END of class FunctionBasedPathConversionTool

}; //namespace
//=============================================================================
//=============================================================================

#endif // _FunctionBasedPathConversionTool_h__


