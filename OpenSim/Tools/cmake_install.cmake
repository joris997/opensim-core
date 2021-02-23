# Install script for directory: /home/none/opensim-core-joris/OpenSim/Tools

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/none/opensim-core-joris-install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosimTools.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosimTools.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosimTools.so"
         RPATH "$ORIGIN/")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/none/opensim-core-joris/libosimTools.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosimTools.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosimTools.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosimTools.so"
         OLD_RPATH "/home/none/opensim-core-joris:/home/none/simbody/lib:"
         NEW_RPATH "$ORIGIN/")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libosimTools.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/OpenSim/Tools" TYPE FILE FILES
    "/home/none/opensim-core-joris/OpenSim/Tools/ActuatorForceTarget.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/ActuatorForceTargetFast.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/AnalyzeTool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/CMC.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/CMCTool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/CMC_Joint.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/CMC_Orientation.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/CMC_Point.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/CMC_Task.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/CMC_TaskSet.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/CorrectionController.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/DynamicsTool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/ForwardTool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/GenericModelMaker.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/IKCoordinateTask.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/IKMarkerTask.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/IKTask.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/IKTaskSet.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/IMUInverseKinematicsTool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/InverseDynamicsTool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/InverseKinematicsTool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/InverseKinematicsToolBase.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/MarkerPair.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/MarkerPairSet.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/MarkerPlacer.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/Measurement.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/MeasurementSet.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/ModelScaler.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/MuscleStateTrackingTask.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/RRATool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/RegisterTypes_osimTools.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/SMC_Joint.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/ScaleTool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/StateTrackingTask.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/Tool.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/TrackingTask.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/VectorFunctionForActuators.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/osimTools.h"
    "/home/none/opensim-core-joris/OpenSim/Tools/osimToolsDLL.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/none/opensim-core-joris/OpenSim/Tools/Test/cmake_install.cmake")

endif()

