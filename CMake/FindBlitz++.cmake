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

IF(Blitzpp_INCLUDE_DIRS)
  FIND_PATH(BLITZ_INCLUDE_DIR blitz/blitz.h  ${Blitzpp_INCLUDE_DIRS})
ELSE(Blitzpp_INCLUDE_DIRS)
  FIND_PATH(BLITZ_INCLUDE_DIR blitz/blitz.h ${BLITZ_HOME}/include $ENV{BLITZ_HOME}/include $ENV{HOME}/include)
  FIND_LIBRARY(BLITZ_LIBRARIES ${Libblitz} ${BLITZ_HOME}/lib $ENV{BLITZ_HOME}/lib $ENV{HOME}/lib)
ENDIF(Blitzpp_INCLUDE_DIRS)

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
