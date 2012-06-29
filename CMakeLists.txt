macro (add_sources)
    file (RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
    foreach (_src ${ARGN})
        if (_relPath)
            list (APPEND SRCS "${_relPath}/${_src}")
        else()
            list (APPEND SRCS "${_src}")
        endif()
    endforeach()
    if (_relPath)
        # propagate SRCS to parent directory
        set (SRCS ${SRCS} PARENT_SCOPE)
    endif()
endmacro()





cmake_minimum_required (VERSION 2.6) 
project (pimc++) 
SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS ,"cc")
SET(CMAKE_CXX_COMPILER,"g++")
SET(CMAKE_CXX_FLAGS "-DH5_USE_16_API -O3 -I/Users/bkclark/Code/MauriceTest/blitz/ -L/Users/bkclark/Code/MauriceTest/blitz/lib  -I/home/bclark/pcts/include/ -I/home/bclark/Install/sprng-2.0/src/ -Wno-write-strings  -DH5_USE_16_API -I/Users/bkclark/Install/hdf5-1.8.7-mac-intel-x86_64-shared/include/  -I/home/bclark/Install/sprng-2.0/src/ -Wno-write-strings -L/home/bclark/pcts/lib/   -L/Users/bkclark/Install/hdf5-1.8.7-mac-intel-x86_64-shared/lib/ -lhdf5 -lsprng -L/usr/lib -L/home/bclark/lib/  -lz  -L/Applications/sage/local/lib/ -I/usr/include/atlas/  -lblas -llapack -lgsl -L/opt/intel/mkl/10.2.6.038/lib/em64t/  -lpthread   ")
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
src/Actions/DavidLongRangeClassYk2.cc																       

src/Moves/OpenStage.cc
src/Actions/DiagonalActionClass.cc
src/SwitchClass.cc						        
src/Moves/PermuteStage.cc
src/Actions/EAMClass.cc																			     
src/Moves/PermuteTableClass.cc
src/Actions/ExternalPotential.cc																	       
src/Moves/PermuteTableOnClass.cc
src/Actions/FreeNodalActionClass.cc																	         
src/Moves/PreSampling.cc
src/Actions/KineticClass.cc																			    
src/Moves/RandomPermClass.cc 																		      
src/Moves/ReadPath.cc																		        
src/Moves/RefSliceMove.cc
src/Actions/LongRangeClass.cc																				   
src/Moves/SPS.cc
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
src/MoleculeHelper.cc																			
src/Observables/TimeLindenman.cc
src/Moves/AVBMove.cc													
src/Moves/BisectionBlock.cc																	
src/Observables/Weight.cc
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
src/WriteData.cc)
message(STATUS "${SRCS}")
add_executable (pimc++ ${SRCS})

