project (pimc++ CXX C Fortran)

cmake_minimum_required (VERSION 2.6)

#-------------------------------------------------------------------
# Options
#-------------------------------------------------------------------

SET(NDIM 3)
SET(OPEN_LOOPS FALSE)
SET(USE_MPI TRUE)
SET(QMC_BUILD_STATIC FALSE)

#-------------------------------------------------------------------
# Compilers
#-------------------------------------------------------------------
SET(C_FLAGS "-DH5_USE_16_API -O3 -msse3 -w -DNDIM=${NDIM}")
IF(OPEN_LOOPS)
  SET(C_FLAGS "${C_FLAGS} -DOPEN_LOOPS=1")
ELSE(OPEN_LOOPS)
  SET(C_FLAGS "${C_FLAGS} -DOPEN_LOOPS=0")
ENDIF(OPEN_LOOPS)
IF(USE_MPI)
  SET(C_FLAGS "${C_FLAGS} -DUSE_MPI -fopenmp")
ENDIF(USE_MPI)


SET(CMAKE_CXX_COMPILER,"CC")
SET(CMAKE_CXX_FLAGS "${C_FLAGS}")
SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS ,"cc")

SET(CMAKE_C_COMPILER,"cc")
SET(CMAKE_C_FLAGS "-std=c99 ${C_FLAGS}")
SET(CMAKE_C_SOURCE_FILE_EXTENSIONS ,"c")

SET(CMAKE_Fortran_COMPILER,"ftn")
SET(CMAKE_Fortran_FLAGS "-funroll-loops -O3")
SET(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS ,"f")

IF(QMC_BUILD_STATIC)
  # Enable static linking
  set(BUILD_STATIC_LIBS TRUE)
  set(LINK_SEARCH_END_STATIC TRUE)
  set(CMAKE_EXE_LINKER_FLAGS "-static -Wl,-Bstatic")
  set(MPI_LINK_FLAGS "-Wl,-Bstatic")
  set(LINK_SEARCH_END_STATIC TRUE)
  set(HDF5_USE_STATIC_LIBRARIES TRUE)

  # Disable shared linking
  set(BUILD_SHARED_LIBS FALSE)
  set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
  set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
  set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
ENDIF(QMC_BUILD_STATIC)

#-------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------
SET(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "Single output directory for building all libraries")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Single output directory for building all executables")
SET(CMAKE_INSTALL_PREFIX "$ENV{HOME}" CACHE PATH "Default install to home directory")
SET(PROJECT_CMAKE "${CMAKE_SOURCE_DIR}/CMake")
SET(BUILD_DIR "${CMAKE_SOURCE_DIR}/build")

#-------------------------------------------------------------------
#  check MPI
#-------------------------------------------------------------------
IF(USE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
    MESSAGE(STATUS "MPI_COMPILE_FLAGS=${MPI_COMPILE_FLAGS}")
    MESSAGE(STATUS "MPI_LINK_FLAGS=${MPI_LINK_FLAGS}")
    MESSAGE(STATUS "MPI_INCLUDE_PATH=${MPI_INCLUDE_PATH}")
    MESSAGE(STATUS "MPI_LIBRARIES=${MPI_LIBRARIES}")
    SET(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
    LINK_LIBRARIES(${MPI_LIBRARIES})
  ELSE(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
    SET(TRIAL_LIBRARY_PATHS
      $ENV{MPI_HOME}/lib
      /usr/lib
      /usr/local/lib
      /opt/lib
      /sw/lib
    )
    SET(TRIAL_INCLUDE_PATHS
      $ENV{MPI_HOME}/include
      /usr/include
      /opt/include
      /usr/local/include
      /sw/include
    )
    FIND_LIBRARY(MPI_LIBRARIES mpi ${TRIAL_LIBRARY_PATHS})
    FIND_PATH(MPI_INCLUDE_PATH mpi.h ${TRIAL_INCLUDE_PATHS})
    IF(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
      MESSAGE(STATUS "MPI_COMPILE_FLAGS=${MPI_COMPILE_FLAGS}")
      MESSAGE(STATUS "MPI_LINK_FLAGS=${MPI_LINK_FLAGS}")
      MESSAGE(STATUS "MPI_INCLUDE_PATH=${MPI_INCLUDE_PATH}")
      MESSAGE(STATUS "MPI_LIBRARIES=${MPI_LIBRARIES}")
      SET(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
      SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
      INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
      LINK_LIBRARIES(${MPI_LIBRARIES})
    ELSE(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
      MESSAGE(STATUS "MPI not found!")
    ENDIF(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
  ENDIF(MPI_INCLUDE_PATH AND MPI_LIBRARIES)
ENDIF(USE_MPI)

#-------------------------------------------------------------------
# check blitz++
#-------------------------------------------------------------------
INCLUDE(${PROJECT_CMAKE}/FindBlitz++.cmake)
IF(BLITZ_FOUND)
  SET(HAVE_BLITZ 1)
  INCLUDE_DIRECTORIES(${BLITZ_INCLUDE_DIR})
  IF(QMC_BUILD_STATIC)
    LINK_LIBRARIES(${BLITZ_LIBRARIES})
  ENDIF(QMC_BUILD_STATIC)
ELSE(BLITZ_FOUND)
  MESSAGE(FATAL_ERROR "Require Blitz++. Set BLITZ_HOME")
ENDIF(BLITZ_FOUND)

#-------------------------------------------------------------------
# check fftw
#-------------------------------------------------------------------
INCLUDE(${PROJECT_CMAKE}/FindFFTW.cmake)
IF(FFTW_FOUND)
  SET(HAVE_FFTW 1)
  INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})
  LINK_LIBRARIES(${FFTW_LIBRARIES})
ELSE(FFTW_FOUND)
  MESSAGE(FATAL_ERROR "Require FFTW. Set FFTW_HOME")
ENDIF(FFTW_FOUND)

#-------------------------------------------------------------------
# check sprng
#-------------------------------------------------------------------
INCLUDE(${PROJECT_CMAKE}/FindSPRNG.cmake)
IF(SPRNG_FOUND)
  SET(HAVE_SPRNG 1)
  INCLUDE_DIRECTORIES(${SPRNG_INCLUDE_DIR})
  IF(QMC_BUILD_STATIC)
    LINK_LIBRARIES(${SPRNG_LIBRARY})
  ELSE(QMC_BUILD_STATIC)
    LINK_LIBRARIES(${SPRNG_LIBRARY} -lgmp)
  ENDIF(QMC_BUILD_STATIC)
ELSE(SPRNG_FOUND)
  MESSAGE(FATAL_ERROR "Require SPRNG. Set SPRNG_HOME")
ENDIF(SPRNG_FOUND)

#-------------------------------------------------------------------
# check gsl
#-------------------------------------------------------------------
INCLUDE(${PROJECT_CMAKE}/FindGSL.cmake)
IF(GSL_FOUND)
  SET(HAVE_LIBGSL 1)
  INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})
  LINK_LIBRARIES(${GSL_LIBRARY})
ELSE(GSL_FOUND)
  MESSAGE(FATAL_ERROR "Require GSL. Set GSL_HOME")
ENDIF(GSL_FOUND)

#-------------------------------------------------------------------
# check hdf5
#-------------------------------------------------------------------
INCLUDE(${PROJECT_CMAKE}/FindHDF5.cmake)
IF(HDF5_FOUND)
  SET(HAVE_LIBHDF5 1)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
  LINK_LIBRARIES(${HDF5_LIBRARIES})
ELSE(HDF5_FOUND)
  MESSAGE(FATAL_ERROR "Require hdf5 1.6.4 or higher. Set HDF5_HOME")
ENDIF(HDF5_FOUND)

#-------------------------------------------------------------------
# check blas/lapack/atlas/mkl
#-------------------------------------------------------------------
INCLUDE(${PROJECT_CMAKE}/FindLapack.cmake)
IF(LAPACK_FOUND AND BLAS_FOUND)
  SET(HAVE_LAPACK 1)
  SET(HAVE_BLAS 1)
ELSE(LAPACK_FOUND AND BLAS_FOUND)
  MESSAGE(FATAL_ERROR "Require lapack/blas/atlas/mkl. Set MKL_HOME etc.")
ENDIF(LAPACK_FOUND AND BLAS_FOUND)

#-------------------------------------------------------------------
# extra libraries
#-------------------------------------------------------------------
#LINK_LIBRARIES("-lz -lpthread -lm -lmpi -mp")
#LINK_LIBRARIES("-lz -lpthread -lm -lsz")
LINK_LIBRARIES("-lz -lpthread -lm")

#-------------------------------------------------------------------
# add sources
#-------------------------------------------------------------------
macro (add_sources)
    file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND SRCS "${_relPath}/${_src}")
        else(_relPath)
            list (APPEND SRCS "${_src}")
        endif(_relPath)
    endforeach(_src ${ARGN})
    if (_relPath)
        # propagate SRCS to parent directory
        set (SRCS ${SRCS} PARENT_SCOPE)
    endif(_relPath)
endmacro (add_sources)

# sources
add_sources(src/Actions/ActionBase.cc
src/Moves/MoveBase.cc
src/Actions/ActionsClass.cc
src/Actions/Sal.cc
src/Moves/MoveUtils.cc
src/Actions/BlendActions.cc
src/Moves/MultiStage.cc
src/Actions/CummingsWaterPotential.cc
src/Moves/NoPermuteClass.cc
src/Actions/DavidLongRangeClass.cc
src/Moves/OpenEndMove.cc
src/Actions/DavidLongRangeClassYk.cc
src/Moves/OpenStage.cc
src/Actions/DiagonalActionClass.cc
src/Moves/PermuteStage.cc
src/Actions/EAMClass.cc
src/Moves/PermuteTableClass.cc
src/Actions/HarmonicPotential.cc
src/Moves/PermuteTableOnClass.cc
src/Actions/FreeNodalActionClass.cc
src/Actions/ParametrizedFreeNodalActionClass.cc
src/Moves/PreSampling.cc
src/Actions/KineticClass.cc
src/Moves/RandomPermClass.cc
src/Moves/ReadPath.cc
src/Moves/RefSliceMove.cc
src/Actions/LongRangeClass.cc
src/Moves/SPS.cc
src/SwitchClass.cc
src/Actions/LongRangeCoulombClass.cc
src/Moves/SpaceWarp.cc
src/Actions/LongRangePotClass.cc
src/Moves/StageClass.cc
src/Actions/LongRangeRPAClass.cc
src/Moves/StructureRejectStage.cc
src/Actions/MoleculeInteractionsClass.cc
src/Moves/SwapMove.cc
src/Actions/Mu.cc
src/Moves/TablePermuteStage.cc
src/Actions/NodalActionClass.cc
src/Actions/OpenLoopImportance.cc
src/Moves/WormClose.cc
src/Actions/QMCSamplingClass.cc
src/Moves/WormGrow.cc
src/Actions/ReadAction.cc
src/Moves/WormPermuteStage.cc
src/Actions/SHONodalActionClass.cc
src/Moves/WormRemove.cc
src/Actions/ST2WaterClass.cc
src/Moves/WormStage.cc
src/Actions/ShortRangeApproximateClass.cc
src/Observables/Angular.cc
src/Actions/ShortRangeClass.cc
src/Observables/Coupling.cc
src/Actions/ShortRangeOnClass.cc
src/Observables/DistanceToHead.cc
src/Actions/WaterClass.cc
src/Actions/WaterFast.cc
src/Actions/dVec.cc
src/Actions/ShortRangeOn_diagonal_Class.cc
src/Observables/DynamicStructureFactor.cc
src/Actions/ShortRangeOn_diagonal_displace_Class.cc
src/Observables/Forces.cc
src/Actions/ShortRangePotClass.cc
src/Observables/HBond.cc
src/Actions/ShortRangePrimitive.cc
src/Observables/Hexatic.cc
src/Actions/StructureReject.cc
src/Observables/JosephsonPathDump.cc
src/Actions/TIP5PWaterClass.cc
src/Observables/ObservableBase.cc
src/Actions/Tether.cc
src/Observables/ObservableCorrelation.cc
src/Observables/ObservableRefCorrelation.cc
src/Observables/ObservableDiffusion.cc
src/Observables/ObservableEnergy.cc
src/Common.cc
src/Observables/ObservableVar.cc
src/Communication/Communication.cc
src/Observables/OpenOrientation.cc
src/EventClass.cc
src/Observables/PairCorrelationReweighting.cc
src/GridClass.cc
src/Observables/ParticleAverageLoc.cc
src/IO/IO.cc
src/Observables/PathDump.cc
src/IO/IOASCII.cc
src/Observables/PermutationCount.cc
src/IO/IOHDF5.cc
src/Observables/Phik.cc
src/IO/IOVarHDF5.cc
src/Observables/PlaneDensity.cc
src/IO/InputFile.cc
src/Observables/Pressure.cc
src/IO/InputOutput.cc
src/Observables/SpecificHeat.cc
src/IO/InputOutputASCII.cc
src/Observables/SpecificHeatA.cc
src/IO/InputOutputHDF5.cc
src/Observables/StructureFactor.cc
src/Observables/RefStructureFactor.cc
src/InitPaths.cc
src/Observables/SuperfluiDrop.cc
src/LoopClass.cc
src/Observables/SuperfluidFraction.cc
src/Observables/Sign.cc
src/Main.cc
src/Observables/SuperfluidFractionPerLayer.cc
src/MatrixOps/MatrixOps.cc
src/Observables/Time.cc
src/MirroredClass.cc
src/Observables/TimeHexatic.cc
src/Observables/Centroid.cc
src/Observables/NodalModelTime.cc
src/Moves/CentroidMove.cc
src/MoleculeHelper.cc
src/Observables/TimeLindenman.cc
src/Moves/AVBMove.cc
src/Moves/BisectionBlock.cc
src/Observables/Weight.cc
src/Observables/ImportanceWeight.cc
src/Moves/BisectionJosephsonStage.cc
src/Observables/WindingNumber.cc
src/Moves/BisectionSphereBlock.cc
src/PIMCClass.cc
src/Moves/BisectionSphereStage.cc
src/PairAction/DavidPAClass.cc
src/Moves/BisectionStage.cc
src/PairAction/PAclassicalFit.cc
src/Moves/CenterDropletMove.cc
src/PairAction/PAzeroFit.cc
src/Moves/CenterofMassMove.cc
src/PairAction/Potential.cc
src/Moves/ClusterMove.cc
src/ParseCommand.cc
src/Moves/CorrelatedBisectionBlock.cc
src/PathClass.cc
src/Moves/CoupledPermuteStage.cc
src/PathDataClass.cc
src/Moves/CouplingMove.cc
src/QMCWrapper.cc
src/Moves/CouplingStage.cc
src/SpeciesClass.cc
src/Moves/DisplaceMove.cc
src/Moves/NodalModelMove.cc
src/Moves/DisplaceMoveFast.cc
src/Moves/EmptyStage.cc
src/Splines/CubicSpline.cc
src/Splines/BicubicSpline.cc
src/Splines/QuinticSpline.cc
src/Moves/EndStage.cc
src/Moves/IonMoveManager.cc
src/Moves/IonMoveStage.cc
src/Splines/Grid.cc
src/Moves/LeviFlightStage.cc
src/Moves/MetaMoves.cc
src/Moves/MoleculeBias.cc
src/Moves/MoleculeMove.cc
src/Moves/MoleculeMoveBase.cc
src/WormHelper.cc
src/Moves/MoleculeMoveManager.cc
src/WriteData.cc
src/Ewald/OptimizedBreakup.cc
src/Integration/GKIntegration.cc
)

message(STATUS "pimc++:${SRCS}")
add_executable (pimc++ ${SRCS})
SET_TARGET_PROPERTIES("pimc++" PROPERTIES LINKER_LANGUAGE CXX)

# sources
SET(SRCS "")
add_sources(src/Ewald/ewald.cc
src/Ewald/OptimizedBreakup.cc
src/MatrixOps/MatrixOps.cc
src/Integration/GKIntegration.cc
src/Integration/HermiteQuad.cc
src/Integration/Integrate.cc
src/IO/IO.cc
src/IO/InputOutput.cc
src/IO/IOASCII.cc
src/IO/IOVarHDF5.cc
src/IO/InputOutputASCII.cc
src/IO/IOHDF5.cc
src/IO/InputFile.cc
src/IO/InputOutputHDF5.cc
src/Splines/BicubicSpline.cc
src/Splines/CubicSpline.cc
src/Splines/Grid.cc
src/Splines/QuinticSpline.cc
)

message(STATUS "ewald:${SRCS}")
add_executable (ewald ${SRCS})
SET_TARGET_PROPERTIES("ewald" PROPERTIES LINKER_LANGUAGE CXX)

# sources
SET(SRCS "")
add_sources(src/FreeParticles.cc
src/MatrixOps/MatrixOps.cc
src/IO/IO.cc
src/Communication/Communication.cc
src/IO/InputOutput.cc
src/IO/IOASCII.cc
src/IO/IOVarHDF5.cc
src/IO/InputOutputASCII.cc
src/IO/IOHDF5.cc
src/IO/InputFile.cc
src/IO/InputOutputHDF5.cc
)

message(STATUS "FreeParticles:${SRCS}")
add_executable (FreeParticles ${SRCS})
SET_TARGET_PROPERTIES("FreeParticles" PROPERTIES LINKER_LANGUAGE CXX)

## sources
#SET(SRCS "")
#add_sources(src/Fitting/Ellipsoid.cc
#src/MatrixOps/MatrixOps.cc
#)
#
#message(STATUS "Ellipsoid:${SRCS}")
#add_executable (Ellipsoid ${SRCS})
#SET_TARGET_PROPERTIES("Ellipsoid" PROPERTIES LINKER_LANGUAGE CXX)

## sources
#SET(SRCS "")
#add_sources(src/test.cc
#)
#
#message(STATUS "test:${SRCS}")
#add_executable (test ${SRCS})
#SET_TARGET_PROPERTIES("test" PROPERTIES LINKER_LANGUAGE CXX)

#EXECUTE_PROCESS( WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
#                 COMMAND chown -R brown281:qmc pimc++ FreeParticles ewald
#                 COMMAND chmod -R 0770 pimc++ FreeParticles ewald )
