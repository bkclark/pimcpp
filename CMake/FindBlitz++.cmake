#
# this module look for blitz++ (http://www.oonumerics.org/blitz) support
# it will define the following values
#
# BLITZ_INCLUDE_DIR = where blitz/blitz.h can be found
#
# May want to define this but seldom required
# BLITZ_LIBRARY = where blitz library can be found (reserved)
#
#IF(EXISTS ${PROJECT_CMAKE}/BlitzppConfig.cmake)
#  INCLUDE(${PROJECT_CMAKE}/BlitzppConfig.cmake)
#ENDIF(EXISTS ${PROJECT_CMAKE}/BlitzppConfig.cmake)

SET(Libblitz blitz)
IF(QMC_BUILD_STATIC)
  SET(Libblitz libblitz.a)
ENDIF(QMC_BUILD_STATIC)

SET(TRIAL_LIBRARY_PATHS
  /usr/local/lib
  /sw/lib
  ${CMAKE_SOURCE_DIR}/lib
  $ENV{BLITZ_HOME}/lib
)

SET(TRIAL_INCLUDE_PATHS
  /usr/local/include
  /sw/include
  ${CMAKE_SOURCE_DIR}/include
  $ENV{BLITZ_HOME}/include
)

FIND_PATH(BLITZ_INCLUDE_DIR blitz/blitz.h ${TRIAL_INCLUDE_PATHS})
FIND_LIBRARY(BLITZ_LIBRARIES ${Libblitz} ${TRIAL_LIBRARY_PATHS})

IF(BLITZ_INCLUDE_DIR)
  SET(BLITZ_FOUND 1 CACHE BOOL "Found blitz++ library")
  MESSAGE(STATUS "BLITZ_INCLUDE_DIR=${BLITZ_INCLUDE_DIR}")
  MESSAGE(STATUS "BLITZ_LIBRARIES=${BLITZ_LIBRARIES}")
ELSE(BLITZ_INCLUDE_DIR)
  SET(BLITZ_FOUND 0 CACHE BOOL "Found blitz++ library")
ENDIF(BLITZ_INCLUDE_DIR)

MARK_AS_ADVANCED(
  BLITZ_INCLUDE_DIR
  BLITZ_FOUND
)
