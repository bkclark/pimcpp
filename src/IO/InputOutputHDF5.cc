/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2013  B. Clark, K. Esler, E. Brown   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
// http://code.google.com/p/pimcplusplus/                  //
/////////////////////////////////////////////////////////////

#include "InputOutput.h"
//#include "InputOutputASCII.h"

#define MAX_HDF5_STRING_LENGTH 200

/************************************************************
 *                    ReadInto Functions                    *
 ************************************************************/

bool VarHDF5Class::ReadInto (double &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dimensions.size() == 1 && Dimensions(0)==1);
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, &val);
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<double,1> &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dimensions.size() == 1);
  val.resize(Dimensions(0));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<double,2> &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dimensions.size() == 2);
  val.resize(Dimensions(0), Dimensions(1));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<double,3> &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dimensions.size() == 3);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}
bool VarHDF5Class::ReadInto (blitz::Array<double,4> &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dimensions.size() == 4);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2), Dimensions(3));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_DOUBLE, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}
 
bool VarHDF5Class::ReadInto (int &val)
{
  assert (Type == INT_TYPE);
  assert (Dimensions.size() == 1 && Dimensions(0)==1);
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, &val);
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<int,1> &val)
{
  assert (Type == INT_TYPE);
  assert (Dimensions.size() == 1);
  val.resize(Dimensions(0));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<int,2> &val)
{
  assert (Type == INT_TYPE);
  assert (Dimensions.size() == 2);
  val.resize(Dimensions(0), Dimensions(1));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<int,3> &val)
{
  assert (Type == INT_TYPE);
  assert (Dimensions.size() == 3);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<int,4> &val)
{
  assert (Type == INT_TYPE);
  assert (Dimensions.size() == 4);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2), Dimensions(3));
  herr_t status = H5Dread(DataSetID, H5T_NATIVE_INT, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, val.data());
  return (status == 0);
}

bool VarHDF5Class::ReadInto (string &val)
{
  assert (Type == STRING_TYPE);
  assert (Dimensions.size() == 1 && Dimensions(0)==1);
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  blitz::Array<char,1> charArray(length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  val = charArray.data();
  H5Tclose(type);
  return (status == 0);
}


bool VarHDF5Class::ReadInto (blitz::Array<string,1> &val)
{
  assert (Type == STRING_TYPE);
  assert (Dimensions.size() == 1);
  val.resize(Dimensions(0));
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  blitz::Array<char,2> charArray(Dimensions(0),length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  for (int i=0; i<Dimensions(0);i++)
    val(i) = &(charArray(i,0));
  H5Tclose(type);
  return (status == 0);
}


bool VarHDF5Class::ReadInto (blitz::Array<string,2> &val)
{
  assert (Type == STRING_TYPE);
  assert (Dimensions.size() == 2);
  val.resize(Dimensions(0), Dimensions(1));
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  blitz::Array<char,3> charArray(Dimensions(0),Dimensions(1),length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  for (int i=0; i<Dimensions(0);i++)
    for (int j=0; j<Dimensions(1);j++)
      val(i,j) = &(charArray(i,j,0));
  H5Tclose(type);
  return (status == 0);
}



bool VarHDF5Class::ReadInto (blitz::Array<string,3> &val)
{
  assert (Type == STRING_TYPE);
  assert (Dimensions.size() == 3);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2));
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  blitz::Array<char,4> charArray(Dimensions(0),Dimensions(1),
			  Dimensions(2),length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  for (int i=0; i<Dimensions(0);i++)
    for (int j=0; j<Dimensions(1);j++)
      for (int k=0; k<Dimensions(2);k++)
	val(i,j,k) = &(charArray(i,j,k,0));
  H5Tclose(type);
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<string,4> &val)
{
  assert (Type == STRING_TYPE);
  assert (Dimensions.size() == 4);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2), Dimensions(3));
  hid_t type = H5Dget_type(DataSetID);
  size_t length = H5Tget_size(type);

  blitz::Array<char,5> charArray(Dimensions(0),Dimensions(1),
			  Dimensions(2),Dimensions(3),length);
  herr_t status = H5Dread(DataSetID, type, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, charArray.data());
  for (int i=0; i<Dimensions(0);i++)
    for (int j=0; j<Dimensions(1);j++)
      for (int k=0; k<Dimensions(2);k++)
	for (int l=0; l<Dimensions(3);l++)
	  val(i,j,k,l) = &(charArray(i,j,k,l,0));
  H5Tclose(type);
  return (status == 0);
}


bool VarHDF5Class::ReadInto (bool &val)
{ 
  assert (Type == BOOL_TYPE);
  assert (Dimensions.size() == 1 && Dimensions(0)==1);
  unsigned char rval;
  herr_t status = H5Dread(DataSetID, BoolType, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, &rval);
  val = (rval != 0);
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<bool,1> &val)
{
  assert (Type == BOOL_TYPE);
  assert (Dimensions.size() == 1);
  val.resize(Dimensions(0));
  val = false;
  blitz::Array<unsigned char,1> rval(val.size());
  herr_t status = H5Dread(DataSetID, BoolType, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, rval.data());
  for (int i=0; i<val.size(); i++)
    val(i) = (rval(i) != 0);
  return (status == 0);
}


bool VarHDF5Class::ReadInto (blitz::Array<bool,2> &val)
{
  assert (Type == INT_TYPE);
  assert (Dimensions.size() == 2);
  val.resize(Dimensions(0), Dimensions(1));
  val = false;
  blitz::Array<unsigned char,2> rval(Dimensions(0), Dimensions(1));
  herr_t status = H5Dread(DataSetID, BoolType, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, rval.data());
  for (int i=0; i<Dimensions(0); i++)
    for (int j=0; j<Dimensions(1); j++)
      val(i,j) = (rval(i,j) != 0);
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<bool,3> &val)
{
  assert (Type == INT_TYPE);
  assert (Dimensions.size() == 3);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2));
  blitz::Array<unsigned char,3> rval(Dimensions(0), Dimensions(1), Dimensions(2));
  val = false;
  herr_t status = H5Dread(DataSetID, BoolType, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, rval.data());
  for (int i=0; i<Dimensions(0); i++)
    for (int j=0; j<Dimensions(1); j++)
      for (int k=0; k<Dimensions(2); k++)
	val(i,j,k) = (rval(i,j,k) != 0);
  return (status == 0);
}

bool VarHDF5Class::ReadInto (blitz::Array<bool,4> &val)
{
  assert (Type == BOOL_TYPE);
  assert (Dimensions.size() == 4);
  val.resize(Dimensions(0), Dimensions(1), Dimensions(2), Dimensions(3));
  blitz::Array<unsigned char,4> rval(Dimensions(0), Dimensions(1),
			      Dimensions(2), Dimensions(3));
  val = false;
  herr_t status = H5Dread(DataSetID, BoolType, H5S_ALL,
			  H5S_ALL, H5P_DEFAULT, rval.data());
  for (int i=0; i<Dimensions(0); i++)
    for (int j=0; j<Dimensions(1); j++)
      for (int k=0; k<Dimensions(2); k++)
	for (int l=0; l<Dimensions(3); l++)
	  val(i,j,k,l) = (rval(i,j,k,l) != 0);
  return (status == 0);
}



/************************************************************
 *                    Append Functions                      *
 ************************************************************/

bool VarHDF5Class::Append(double val)
{
  const int RANK = 1;

  if ((Type != DOUBLE_TYPE) || (Dim != 1)) {
    cerr << "Trying to append a double to a non blitz::Array<double,1>.  Exitting.\n";
    return false;
  }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;
  hsize_t maxdim[1];
  maxdim[0] = H5S_UNLIMITED;
  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);
  hsize_t size[1];
  hsize_t offset[1];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  H5Dwrite (DataSetID, H5T_NATIVE_DOUBLE, memSpace, DataSpaceID,
  	    H5P_DEFAULT, &val);
  H5Sclose (memSpace);
  return (true);
}



bool VarHDF5Class::Append(blitz::Array<double,1> &val)
{
  const int RANK = 2;

  if ((Type != DOUBLE_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<double,1> to a non blitz::Array<double,2>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<double,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, H5T_NATIVE_DOUBLE, memSpace, DataSpaceID,
  	    H5P_DEFAULT, val.data());
  H5Sclose (memSpace);
  return (true);
}



bool VarHDF5Class::Append(blitz::Array<double,2> &val)
{
  const int RANK = 3;

  if ((Type != DOUBLE_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<double,2> to a non blitz::Array<double,3>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<double,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, H5T_NATIVE_DOUBLE, memSpace, DataSpaceID,
  	    H5P_DEFAULT, val.data());
  H5Sclose (memSpace);
  return (true);
}


bool VarHDF5Class::Append(blitz::Array<double,3> &val)
{
  const int RANK = 4;

  if ((Type != DOUBLE_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<double,3> to a non blitz::Array<double,4>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<double,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, H5T_NATIVE_DOUBLE, memSpace, DataSpaceID,
  	    H5P_DEFAULT, val.data());
  H5Sclose (memSpace);
  return (true);
}





bool VarHDF5Class::Append(int val)
{
  const int RANK = 1;

  if ((Type != INT_TYPE) || (Dim != 1)) {
    cerr << "Trying to append a int to a non blitz::Array<int,1>.  Exitting.\n";
    return false;
  }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;
  hsize_t maxdim[1];
  maxdim[0] = H5S_UNLIMITED;
  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);
  hsize_t size[1];
  hsize_t offset[1];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  H5Dwrite (DataSetID, H5T_NATIVE_INT, memSpace, DataSpaceID,
  	    H5P_DEFAULT, &val);
  H5Sclose (memSpace);
  return (true);
}



bool VarHDF5Class::Append(blitz::Array<int,1> &val)
{
  const int RANK = 2;

  if ((Type != INT_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<int,1> to a non blitz::Array<int,2>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<int,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, H5T_NATIVE_INT, memSpace, DataSpaceID,
  	    H5P_DEFAULT, val.data());
  H5Sclose (memSpace);
  return (true);
}



bool VarHDF5Class::Append(blitz::Array<int,2> &val)
{
  const int RANK = 3;

  if ((Type != INT_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<int,2> to a non blitz::Array<int,3>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<int,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, H5T_NATIVE_INT, memSpace, DataSpaceID,
  	    H5P_DEFAULT, val.data());
  H5Sclose (memSpace);
  return (true);
}



bool VarHDF5Class::Append(blitz::Array<int,3> &val)
{
  const int RANK = 4;

  if ((Type != INT_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<int,3> to a non blitz::Array<int,4>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<int,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, H5T_NATIVE_INT, memSpace, DataSpaceID,
  	    H5P_DEFAULT, val.data());
  H5Sclose (memSpace);
  return (true);
}



bool VarHDF5Class::Append(string val)
{
  const int RANK = 1;

  if ((Type != STRING_TYPE) || (Dim != 1)) {
    cerr << "Trying to append a string to a non blitz::Array<string,1>.  Exitting.\n";
    return false;
  }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;
  hsize_t maxdim[1];
  maxdim[0] = H5S_UNLIMITED;
  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);
  hsize_t size[1];
  hsize_t offset[1];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;

  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, MAX_HDF5_STRING_LENGTH+1);


  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  H5Dwrite (DataSetID, strType, memSpace, DataSpaceID,
  	    H5P_DEFAULT, val.c_str());
  H5Sclose (memSpace);
  return (true);
}



bool VarHDF5Class::Append(blitz::Array<string,1> &strs)
{
  const int RANK = 2;
  if ((Type != STRING_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<string,1> to a non blitz::Array<string,2>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=strs.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<string,2>.\n";
      return false;
    }

  // Copy strings over into a character array
  blitz::Array<char,2> charArray(strs.extent(0), 
			  MAX_HDF5_STRING_LENGTH+1);

  for (int i=0; i<strs.extent(0); i++)
    {
      for (int x=0; 
	   (x<strs(i).length()) && (x < MAX_HDF5_STRING_LENGTH); 
	   x++) 
	charArray(i,x) = (strs(i))[x];
      int x = strs(i).length();
      if (x > MAX_HDF5_STRING_LENGTH)
	x = MAX_HDF5_STRING_LENGTH;
      charArray(i,x) = '\0';
    }
  

  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, MAX_HDF5_STRING_LENGTH+1);


  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }


  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, strType, memSpace, DataSpaceID,
  	    H5P_DEFAULT,charArray.data());
  H5Sclose (memSpace);
  return (true);
}

bool VarHDF5Class::Append(blitz::Array<string,2> &strs)
{
  const int RANK = 3;
  if ((Type != STRING_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<string,2> to a non blitz::Array<string,3>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=strs.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<string,2>.\n";
      return false;
    }

  // Copy strings over into a character array
  blitz::Array<char,3> charArray(strs.extent(0),strs.extent(1),
			  MAX_HDF5_STRING_LENGTH+1);




  for (int i=0; i<strs.extent(0); i++)
    for (int j=0; j<strs.extent(1); j++)
      {
	for (int x=0; 
	     (x<strs(i,j).length()) && (x < MAX_HDF5_STRING_LENGTH); 
	     x++) 
	  charArray(i,j,x) = (strs(i,j))[x];
	int x = strs(i,j).length();
	if (x > MAX_HDF5_STRING_LENGTH)
	  x = MAX_HDF5_STRING_LENGTH;
	charArray(i,j,x) = '\0';
      }
  
  
  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, MAX_HDF5_STRING_LENGTH+1);


  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }


  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, strType, memSpace, DataSpaceID,
  	    H5P_DEFAULT,charArray.data());
  H5Sclose (memSpace);
  return (true);
}


bool VarHDF5Class::Append(blitz::Array<string,3> &strs)
{
  const int RANK = 4;
  if ((Type != STRING_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<string,3> to a non blitz::Array<string,4>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=strs.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<string,2>.\n";
      return false;
    }

  // Copy strings over into a character array
  blitz::Array<char,4> charArray(strs.extent(0),strs.extent(1), strs.extent(2),
			  MAX_HDF5_STRING_LENGTH+1);

  for (int i=0; i<strs.extent(0); i++)
    for (int j=0; j<strs.extent(1); j++)
      for (int k=0; k<strs.extent(2); k++)
	{
	  for (int x=0; 
	       (x<strs(i,j,k).length()) && (x < MAX_HDF5_STRING_LENGTH); 
	       x++) 
	    charArray(i,j,k,x) = (strs(i,j,k))[x];
	  int x = strs(i,j,k).length();
	  if (x > MAX_HDF5_STRING_LENGTH)
	    x = MAX_HDF5_STRING_LENGTH;
	  charArray(i,j,k,x) = '\0';
	}
  
  
  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, MAX_HDF5_STRING_LENGTH+1);

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  H5Dwrite (DataSetID, strType, memSpace, DataSpaceID,
  	    H5P_DEFAULT,charArray.data());
  H5Sclose (memSpace);
  return (true);
}



bool VarHDF5Class::Append(bool val)
{
  const int RANK = 1;

  if ((Type != BOOL_TYPE) || (Dim != 1)) {
    cerr << 
      "Trying to append a bool to a non blitz::Array<bool,1>.  Exitting.\n";
    return false;
  }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;
  hsize_t maxdim[1];
  maxdim[0] = H5S_UNLIMITED;
  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);
  hsize_t size[1];
  hsize_t offset[1];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  unsigned char wval = val ? (unsigned char) 1 : (unsigned char)0;
  H5Dwrite (DataSetID, BoolType, memSpace, DataSpaceID,
  	    H5P_DEFAULT, &wval);
  H5Sclose (memSpace);
  return (true);
}

bool VarHDF5Class::Append(blitz::Array<bool,1> &val)
{
  const int RANK = 2;

  if ((Type != BOOL_TYPE) || (Dim != RANK)) {
    cerr << 
      "Trying to append an Array<bool,1> to a non blitz::Array<int,2>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<int,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  blitz::Array<unsigned char,1> wval(val.extent(0));
  for (int i=0; i<val.extent(0); i++)
    wval(i) = val(i) ? (unsigned char) 1 : (unsigned char) 0;
  H5Dwrite (DataSetID, BoolType, memSpace, DataSpaceID,
  	    H5P_DEFAULT, wval.data());
  H5Sclose (memSpace);
  return (true);
}


bool VarHDF5Class::Append(blitz::Array<bool,2> &val)
{
  const int RANK = 3;

  if ((Type != BOOL_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<bool,2> " 
	 << "to a non blitz::Array<bool,3>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<bool,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  blitz::Array<unsigned char,2> wval(val.extent(0),val.extent(1));
  for (int i=0; i<val.extent(0); i++)
    for (int j=0; j<val.extent(1); j++)
      wval(i,j) = val(i,j) ? (unsigned char) 1 : (unsigned char) 0;
  H5Dwrite (DataSetID, BoolType, memSpace, DataSpaceID,
  	    H5P_DEFAULT, wval.data());
  H5Sclose (memSpace);
  return (true);
}

bool VarHDF5Class::Append(blitz::Array<bool,3> &val)
{
  const int RANK = 4;

  if ((Type != BOOL_TYPE) || (Dim != RANK)) {
    cerr << "Trying to append an Array<bool,3> "
	 << "to a non blitz::Array<bool,4>.  Exitting.\n";
    return false;
  }

  for (int i=1; i<RANK; i++)
    if (Dimensions(i)!=val.extent(i-1)) {
      cerr << "Dimension mismatch in appending to blitz::Array<int,2>.\n";
      return false;
    }

  hid_t memSpace;
  // Extend the dataset
  Dimensions(0)++;

  hsize_t maxdim[RANK], size[RANK];
  hsize_t offset[RANK];
  offset[0] = Dimensions(0)-1;
  size[0] = 1;
  maxdim[0] = H5S_UNLIMITED;
  for (int i=1; i<RANK; i++) {
    offset[i] = 0;
    size[i] = Dimensions(i);
    maxdim[i] = Dimensions(i);
  }

  H5Dextend(DataSetID, Dimensions.data());
  H5Sset_extent_simple(DataSpaceID, RANK, Dimensions.data(), maxdim);

  /// Select which part of the file to write to
  H5Sselect_hyperslab(DataSpaceID, H5S_SELECT_SET, offset, NULL,
		      size, NULL);
  memSpace = H5Screate_simple (RANK, size, NULL);
  /// Write the new data
  blitz::Array<unsigned char,3> wval(val.extent(0),val.extent(1),val.extent(2));
  for (int i=0; i<val.extent(0); i++)
    for (int j=0; j<val.extent(1); j++)
      for (int k=0; k<val.extent(2); k++)
	wval(i,j,k) = val(i,j,k) ? (unsigned char) 1 : (unsigned char) 0;
  H5Dwrite (DataSetID, BoolType, memSpace, DataSpaceID,
  	    H5P_DEFAULT, wval.data());
  H5Sclose (memSpace);
  return (true);
}

/************************************************************
 *                     Helper Functions                     *
 ************************************************************/

/// Strips everything after and including a '.' in the string.
/// Used to remove section numbers.
// void IOTreeHDF5Class::StripName (string str,string &newString,
// 				    int &myInt)
// {
//   int pos = str.find(".");
//   //  assert(pos>0);
//   if (pos<=0){
//     myInt=0;
//     newString=str;
//   }
//   else{
//     newString=str.substr(0,pos);
//     string intString=str.substr(pos+1,str.length()-1);
//     char *endptr;
//     myInt=strtol(intString.c_str(),&endptr,10);
//     assert (*endptr=='\0');
//   }
// }

/// Strips everything after and including a '_' in the string.
/// Used to remove section numbers.
void IOTreeHDF5Class::StripName (string str, string &newString,
				 int &myInt)
{
  int pos = str.length()-1;
  while ((pos>0) && (str[pos]!='_'))
    pos--;
  if (pos<=0){
    myInt=0;
    newString=str;
  }
  else{
    newString=str.substr(0,pos);
    string intString=str.substr(pos+1,str.length()-1);
    char *endptr;
    myInt=strtol(intString.c_str(),&endptr,10);
    assert (*endptr=='\0');
  }
}


/// Takes an integer and returns a string which is in the format .# 
string NumExtension (int num)
{
  string retString = "_";
  char numstr[100];
  snprintf(numstr, 100, "%d", num);
  //int len = strlen(numstr);
  //int numZeros = 8-len;
  //for (int i=0; i<numZeros; i++)
  //  retString += "0";
  retString += numstr;
  return (retString);
}


/************************************************************
 *                      File Functions                      *
 ************************************************************/

bool IOTreeHDF5Class::OpenFile(string fileName,
			       string mySectionName,
			       IOTreeClass *parent)
{
  // Turn off error printing
  H5E_auto_t func;
  void *client_data;
  H5Eget_auto(&func, &client_data);
  H5Eset_auto(NULL, NULL);

  GroupID = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  BoolType = H5Topen(GroupID, "BOOL");
  // If we don't already have this in here (old version of this
  // library created this file), create a new one and stick it in the
  // file. 
  if (BoolType < 0) {
      BoolType = H5Tcreate(H5T_ENUM, sizeof(unsigned char));
      unsigned char val = 0;
      H5Tenum_insert(BoolType, "FALSE", &val);
      val = 1;
      H5Tenum_insert(BoolType, "TRUE", &val);
      H5Tcommit (GroupID, "BOOL", BoolType);
  }
  // And turn it back on;
  H5Eset_auto(func, client_data);


  if (GroupID < 0) {
    cerr << "Cannot open file " << fileName << endl;
    return false;
  }

  IsOpen = true;
  Parent = parent;
  ReadGroup (GroupID, "/", parent);
  cerr << "Before name = mySectionName\n";
  Name = mySectionName;
  cerr << "Before FileName=fileName\n";
  FileName=fileName;
  cerr << "Before return\n";
  return true;
}

bool IOTreeHDF5Class::NewFile(string fileName,string myName,IOTreeClass* parent)
{
  FileName = fileName;
  bool success = true;
  
  hid_t FileID = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, 
			   H5P_DEFAULT, H5P_DEFAULT);
  // Create a bool enum type
  BoolType = H5Tcreate(H5T_ENUM, sizeof(unsigned char));
  unsigned char val = (unsigned char) false;
  H5Tenum_insert(BoolType, "FALSE", &val);
  val = (unsigned char) true;
  H5Tenum_insert(BoolType, "TRUE", &val);
  H5Tcommit (FileID, "BOOL", BoolType);

  IsOpen = (FileID >= 0);
  success = IsOpen;
  if (success)
    {
      GroupID = FileID;
      Parent=parent;
      Name=myName;
    }
  return (success);
}


void IOTreeHDF5Class::CloseFile()
{
  // First, free all the variables in the list
  while (!VarList.empty()) {
    delete(VarList.front());
    VarList.pop_front();
  }
   
  // Now, call all closes recursively and delete all sections
  while (!SectionList.empty()) {
      SectionList.front()->CloseFile();
      delete SectionList.front();
      SectionList.pop_front();
  }

  if (FileName!="") {
    // Release BoolType;
    H5Tclose (BoolType);
    H5Fclose(GroupID);
  }
  else {
    H5Gclose(GroupID);
  }
}    


void IOTreeHDF5Class::FlushFile()
{
  // First, flush myself if I'm a root file.
  if (FileName!="") {
    herr_t err = H5Fflush(GroupID, H5F_SCOPE_GLOBAL);
    assert((int)err >= 0);
  }

  list<IOTreeClass*>::iterator iter = SectionList.begin();
  while (iter != SectionList.end()) {
    (*iter)->FlushFile();
    iter++;
  }
}


/************************************************************
 *                    Section Functions                     *
 ************************************************************/

IOTreeClass* IOTreeHDF5Class::NewSection(string newName)
{
  IOTreeHDF5Class* newSection = new IOTreeHDF5Class();
  newSection->Name=newName;
  newSection->Parent=this;
  newSection->MyNumber=CurrSecNum;
  newSection->BoolType = BoolType;
  SectionList.push_back(newSection);
  
  string numstr = NumExtension(CurrSecNum);
  newName += numstr;

  newSection->GroupID = H5Gcreate(GroupID,
				  newName.c_str(), 0);

  CurrSecNum++;
  return newSection;
}

void IOTreeHDF5Class::IncludeSection(IOTreeClass *newSection)
{
  // Assign number
  newSection->MyNumber = CurrSecNum++;
  // Add to section list
  SectionList.push_back(newSection);
  // Put entry in my group in the HDF5 file
  string nameStr = newSection->Name + NumExtension (newSection->MyNumber);
  hid_t newGroupID = H5Gcreate(GroupID,
			       nameStr.c_str(), 0);
  hsize_t dim[1];
  dim[0] = 1;
  hid_t dataspace_id = H5Screate_simple(1, dim, NULL);
  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, newSection->FileName.length()+1);
  // Create a dataset named "include" which has the filename to include.
  hid_t dataset_id =   H5Dcreate(newGroupID, "include", strType, dataspace_id,
				 H5P_DEFAULT);
  
  herr_t status = H5Dwrite(dataset_id, strType, 
			   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			   newSection->FileName.c_str());
  if (status < 0)
    {
      cerr << "Error including to HDF5 file in IncludeSection.\n";
      exit(1);
    }

  H5Dclose (dataset_id);
  H5Sclose (dataspace_id);
  H5Tclose (strType);
  H5Gclose (newGroupID);
}



/// C-style wrapper for member function iterator.
herr_t HDF5GroupIterator(hid_t group_id, const char *member_name,
			 void *classPtr)
{
  cerr << "Before C group iterator.\n";
  IOTreeHDF5Class &HDF5sec= *((IOTreeHDF5Class *)classPtr);
  HDF5sec.GroupIterator(member_name);
  cerr << "After C group iterator.\n";
  return (0);
}


inline AtomicType H5toAtomic (H5T_class_t type)
{
  switch (type) {
  case H5T_FLOAT:
    return DOUBLE_TYPE;
  case H5T_INTEGER:
    return INT_TYPE;
  case H5T_STRING:
    return STRING_TYPE;
  case H5T_ENUM:
    return BOOL_TYPE;
  default:
    return NOT_ATOMIC;
  }
}



void IOTreeHDF5Class::GroupIterator(string member_name)
{
  cerr << "GroupIterator( " << member_name << ")\n";

  H5G_stat_t statbuf;
  
  H5Gget_objinfo(GroupID, member_name.c_str(), 0, &statbuf);
  

  if (statbuf.type == H5G_GROUP) {
    // Open the group
    hid_t newGroupID;
    newGroupID = H5Gopen (GroupID, member_name.c_str());
    if ((int)newGroupID < 0) {
      cerr << "Error in IOTreeHDF5Class::GroupIterator.\n";
      exit(1);
    }
    // First check to see if we're including another file.
    // Turn off error printing
    H5E_auto_t func;
    void *client_data;
    H5Eget_auto(&func, &client_data);
    H5Eset_auto(NULL, NULL);
    hid_t includeData = H5Dopen(newGroupID, "include");
    // And turn it back on;
    H5Eset_auto(func, client_data);
    if ((int)includeData < 0) {
      // We're not including
      IOTreeHDF5Class *newTree = new IOTreeHDF5Class;
      newTree->GroupID = newGroupID;
      StripName(member_name,newTree->Name,newTree->MyNumber);
      InsertSection(newTree);
      newTree->BoolType = BoolType;
      newTree->ReadGroup (GroupID, member_name, this);
    }
    else {
      // We're including
      hid_t type = H5Dget_type(includeData);
      size_t length = H5Tget_size(type);
      blitz::Array<char, 1> charArray(length+1);
      herr_t status = H5Dread(includeData, type, H5S_ALL,
			      H5S_ALL, H5P_DEFAULT, charArray.data());
      string fileName = charArray.data();
      H5Tclose(type);
      H5Dclose(includeData);
      string myName; int myNum;
      StripName (member_name, myName, myNum);
      IOTreeClass *newTree = ReadTree(fileName, myName, this);
      newTree->MyNumber = myNum;
      InsertSection (newTree);
    }
  }
  else if (statbuf.type == H5G_DATASET) {
    VarHDF5Class *newVar = new VarHDF5Class;
    newVar->BoolType = BoolType;
    newVar->DataSetID = H5Dopen(GroupID, member_name.c_str());
    newVar->Name = member_name;
    hid_t dataTypeID = H5Dget_type (newVar->DataSetID);
    newVar->Type = H5toAtomic (H5Tget_class (dataTypeID));
    H5Tclose (dataTypeID);
    newVar->DataSpaceID = H5Dget_space (newVar->DataSetID);
    newVar->Dim = H5Sget_simple_extent_ndims(newVar->DataSpaceID);
    newVar->Dimensions.resize(newVar->Dim);
    H5Sget_simple_extent_dims(newVar->DataSpaceID, newVar->Dimensions.data(), 
			      NULL);
    // If we have a 1D array with a single element, it's really an atomic element,
    // not an array at all.
    if (newVar->Dim == 1 && newVar->Dimensions(0) == 1)
      newVar->Dim = 0;
    
    VarList.push_back(newVar);
  }
  else if (statbuf.type == H5G_TYPE) {
    if (member_name != "BOOL")
      cerr << "Compound types not yet supported "
	   << "in IOTreeHDF5Class.  Ignoring " 
	   << member_name << endl;
  }
  else
    cerr << " Unable to identify an object " << member_name << endl;;
  cerr << "GroupIterator( " << member_name << ") done\n";
}


/************************************************************
 *                    Printing Functions                    *
 ************************************************************/


void PrintIndent(int num)
{
  for (int counter=0;counter<num*2;counter++){
    cout<<' ';
  }
}



void IOTreeHDF5Class::PrintTree(int indentNum)
{
  PrintIndent(indentNum);
  cout<<"Section: "<<Name<<endl;
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    PrintIndent(indentNum+1);
    cout<<"Variable: "<<(*varIter)->Name<<" "<<endl;
    varIter++;
  }
  list<IOTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    (*secIter)->PrintTree(indentNum+1);
    secIter++;
  }
}

void IOTreeHDF5Class::PrintTree()
{
  PrintTree(0);
}

/// ReadGroup iterates over the members of it's group, creating
/// VarHDF5Class objects and new IOTreeHDF5Class objects as it
/// goes, calling itself recursively as necessary to traverse all the
/// subobjects below itself.
void IOTreeHDF5Class::ReadGroup(hid_t parentGroupID,
				string name,
				IOTreeClass *parent)
{
  Parent = parent;
  int listLoc;
  StripName(name,Name,listLoc);  
  
  H5Giterate (parentGroupID, name.c_str(), (int *)NULL, HDF5GroupIterator,
	      this);
  if (SectionList.empty())
    CurrSecNum = 0;
  else
    CurrSecNum = (SectionList.back())->MyNumber+1;
  cerr << "After ReadGroup " << name << endl;
}





/************************************************************
 *                    Output Functions                      *
 ************************************************************/ 

void IOTreeHDF5Class::WriteVar(string name, double T)
{
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = DOUBLE_TYPE;
  newVar->Dim = 0;

  /// Add it to the list
  VarList.push_back(newVar);

  hsize_t dim[1];
  dim[0] = 1;
  // Create the dataspace
  newVar->DataSpaceID = H5Screate_simple(1, dim, NULL);
 
  // Create the dataset
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_DOUBLE, newVar->DataSpaceID,
				H5P_DEFAULT);
  // Write the dataset to the file.
  herr_t status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_DOUBLE, 
			   H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
  if (status < 0)
    cerr << "Error writing double to HDF5 file in WriteVar.\n";
}


void IOTreeHDF5Class::WriteVar(string name, blitz::Array<double,1> &v)
{
  const int RANK=1;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = DOUBLE_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_DOUBLE, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_DOUBLE, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  if (status < 0)
    cerr << "Error writing double to HDF5 file in WriteVar.\n";

  // Free DSprops
  H5Pclose (DSprops);
}


void IOTreeHDF5Class::WriteVar(string name, blitz::Array<double,2> &v)
{
  const int RANK=2;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = DOUBLE_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_DOUBLE, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_DOUBLE, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());

  if (status < 0)
    cerr << "Error writing double to HDF5 file in WriteVar.\n";

  // Free DSprops
  H5Pclose (DSprops);
}



void IOTreeHDF5Class::WriteVar(string name, blitz::Array<double,3> &v)
{
  const int RANK=3;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = DOUBLE_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_DOUBLE, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_DOUBLE, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());

  if (status < 0)
    cerr << "Error writing double to HDF5 file in WriteVar.\n";

  // Free DSprops
  H5Pclose (DSprops);
}



void IOTreeHDF5Class::WriteVar(string name, blitz::Array<double,4> &v)
{
  const int RANK=4;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = DOUBLE_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_DOUBLE, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_DOUBLE, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());

  // Free DSprops
  H5Pclose (DSprops);

  if (status < 0)
    cerr << "Error writing double to HDF5 file in WriteVar.\n";
}




void IOTreeHDF5Class::WriteVar(string name, int T)
{
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = INT_TYPE;
  newVar->Dim = 0;

  /// Add it to the list
  VarList.push_back(newVar);

  hsize_t dim[1];
  dim[0] = 1;
  // Create the dataspace
  newVar->DataSpaceID = H5Screate_simple(1, dim, NULL);
 
  // Create the dataset
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_INT, newVar->DataSpaceID,
				H5P_DEFAULT);
  // Write the dataset to the file.
  herr_t status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_INT, 
			   H5S_ALL, H5S_ALL, H5P_DEFAULT, &T);
  if (status < 0)
    cerr << "Error writing int to HDF5 file in WriteVar.\n";
}


void IOTreeHDF5Class::WriteVar(string name, blitz::Array<int,1> &v)
{
  const int RANK=1;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = INT_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_INT, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_INT, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());

  if (status < 0)
    cerr << "Error writing int to HDF5 file in WriteVar.\n";

  // Free DSprops
  H5Pclose (DSprops);
}


void IOTreeHDF5Class::WriteVar(string name, blitz::Array<int,2> &v)
{
  const int RANK=2;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = INT_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_INT, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_INT, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());

  if (status < 0)
    cerr << "Error writing int to HDF5 file in WriteVar.\n";

  // Free DSprops
  H5Pclose (DSprops);
}



void IOTreeHDF5Class::WriteVar(string name, blitz::Array<int,3> &v)
{
  const int RANK=3;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = INT_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_INT, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_INT, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());

  // Free DSprops
  H5Pclose (DSprops);

  if (status < 0)
    cerr << "Error writing int to HDF5 file in WriteVar.\n";
}



void IOTreeHDF5Class::WriteVar(string name, blitz::Array<int,4> &v)
{
  const int RANK=4;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = INT_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				H5T_NATIVE_INT, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, H5T_NATIVE_INT, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());

  // Free DSprops
  H5Pclose (DSprops);

  if (status < 0)
    cerr << "Error writing int to HDF5 file in WriteVar.\n";
}








void IOTreeHDF5Class::WriteVar(string name,string str)
{
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = STRING_TYPE;
  newVar->Dim = 0;

  /// Add it to the list
  VarList.push_back(newVar);

  hsize_t dim[1];
  dim[0] = 1;
  // Create the dataspace
  newVar->DataSpaceID = H5Screate_simple(1, dim, NULL);
  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, str.length()+1);

  // Create the dataset
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				strType, newVar->DataSpaceID,
				H5P_DEFAULT);
  // Write the dataset to the file.
  herr_t status = H5Dwrite(newVar->DataSetID, strType, 
			   H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			   str.c_str());
  H5Tclose(strType);
  if (status < 0)
    cerr << "Error writing string to HDF5 file in WriteVar.\n";
}




void IOTreeHDF5Class::WriteVar(string name, blitz::Array<string,1> &strs)
{
  const int RANK=1;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = STRING_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = strs.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Copy strings over into a character array
  blitz::Array<char,2> charArray(strs.extent(0), 
			  MAX_HDF5_STRING_LENGTH+1);
  for (int i=0; i<strs.extent(0); i++)
    {
      for (int x=0; 
	   (x<strs(i).length()) && (x < MAX_HDF5_STRING_LENGTH); 
	   x++) 
	charArray(i,x) = (strs(i))[x];
      int x = strs(i).length();
      if (x > MAX_HDF5_STRING_LENGTH)
	x = MAX_HDF5_STRING_LENGTH;
      charArray(i,x) = '\0';
    }
  
  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = strs.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = strs.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, MAX_HDF5_STRING_LENGTH+1);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=strs.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				strType, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, strType, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, charArray.data());
  if (status < 0)
    cerr << "Error writing string to HDF5 file in WriteVar.\n";
  H5Tclose(strType);

  // Free DSprops
  H5Pclose (DSprops);
}





void IOTreeHDF5Class::WriteVar(string name, blitz::Array<string,2> &strs)
{
  const int RANK=2;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = STRING_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = strs.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Copy strings over into a character array
  blitz::Array<char,3> charArray(strs.extent(0), strs.extent(1),
			       MAX_HDF5_STRING_LENGTH+1);
  for (int i=0; i<strs.extent(0); i++)
    for (int j=0; j<strs.extent(1); j++)
    {
      for (int x=0; 
	   (x<strs(i,j).length()) && (x < MAX_HDF5_STRING_LENGTH); 
	   x++) 
	charArray(i,j,x) = (strs(i,j))[x];
      int x = strs(i,j).length();
      if (x > MAX_HDF5_STRING_LENGTH)
	x = MAX_HDF5_STRING_LENGTH;
      charArray(i,j,x) = '\0';
    }
  
  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = strs.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = strs.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, MAX_HDF5_STRING_LENGTH+1);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=strs.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				strType, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, strType, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, charArray.data());
  if (status < 0)
    cerr << "Error writing string to HDF5 file in WriteVar.\n";
  H5Tclose(strType);

  // Free DSprops
  H5Pclose (DSprops);
}




void IOTreeHDF5Class::WriteVar(string name, blitz::Array<string,3> &strs)
{
  const int RANK=3;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = STRING_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = strs.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Copy strings over into a character array
  blitz::Array<char,4> charArray(strs.extent(0), strs.extent(1), strs.extent(2),
			  MAX_HDF5_STRING_LENGTH+1);
  for (int i=0; i<strs.extent(0); i++)
    for (int j=0; j<strs.extent(1); j++)
      for (int k=0; k<strs.extent(2); k++)
	{
	  for (int x=0; 
	       (x<strs(i,j,k).length()) && (x < MAX_HDF5_STRING_LENGTH); 
	       x++) 
	    charArray(i,j,k,x) = (strs(i,j,k))[x];
	  int x = strs(i,j,k).length();
	  if (x > MAX_HDF5_STRING_LENGTH)
	    x = MAX_HDF5_STRING_LENGTH;
	  charArray(i,j,k,x) = '\0';
	}
  
  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = strs.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = strs.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, MAX_HDF5_STRING_LENGTH+1);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=strs.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				strType, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, strType, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, charArray.data());
  if (status < 0)
    cerr << "Error writing string to HDF5 file in WriteVar.\n";
  H5Tclose(strType);

  // Free DSprops
  H5Pclose (DSprops);
}



void IOTreeHDF5Class::WriteVar(string name, blitz::Array<string,4> &strs)
{
  const int RANK=4;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->Name = name;
  newVar->Type = STRING_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = strs.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Copy strings over into a character array
  blitz::Array<char,5> charArray(strs.extent(0), strs.extent(1), strs.extent(2),
			  strs.extent(3), MAX_HDF5_STRING_LENGTH+1);
  for (int i=0; i<strs.extent(0); i++)
    for (int j=0; j<strs.extent(1); j++)
      for (int k=0; k<strs.extent(2); k++)
	for (int l=0; l<strs.extent(3); l++) {
	  for (int x=0; 
	       (x<strs(i,j,k,l).length()) && (x < MAX_HDF5_STRING_LENGTH); 
	       x++) 
	    charArray(i,j,k,l,x) = (strs(i,j,k,l))[x];
	  int x = strs(i,j,k,l).length();
	  if (x > MAX_HDF5_STRING_LENGTH)
	    x = MAX_HDF5_STRING_LENGTH;
	  charArray(i,j,k,l,x) = '\0';
	}
  
  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = strs.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = strs.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  hid_t strType = H5Tcopy (H5T_C_S1);
  H5Tset_size (strType, MAX_HDF5_STRING_LENGTH+1);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=strs.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				strType, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  status = H5Dwrite(newVar->DataSetID, strType, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, charArray.data());
  if (status < 0)
    cerr << "Error writing string to HDF5 file in WriteVar.\n";
  H5Tclose(strType);

  // Free DSprops
  H5Pclose (DSprops);
}





void IOTreeHDF5Class::WriteVar(string name, bool T)
{
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->BoolType = BoolType;
  newVar->Name = name;
  newVar->Type = BOOL_TYPE;
  newVar->Dim = 0;

  /// Add it to the list
  VarList.push_back(newVar);

  hsize_t dim[1];
  dim[0] = 1;
  // Create the dataspace
  newVar->DataSpaceID = H5Screate_simple(1, dim, NULL);
 
  // Create the dataset
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				BoolType, newVar->DataSpaceID,
				H5P_DEFAULT);
  // Write the dataset to the file.
  unsigned char wval = T ? (unsigned char)1 : (unsigned char)0;
  herr_t status = H5Dwrite(newVar->DataSetID, BoolType, 
			   H5S_ALL, H5S_ALL, H5P_DEFAULT, &wval);
  if (status < 0)
    cerr << "Error writing bool to HDF5 file in WriteVar.\n";
}


void IOTreeHDF5Class::WriteVar(string name, blitz::Array<bool,1> &v)
{
  const int RANK=1;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->BoolType = BoolType;
  newVar->Name = name;
  newVar->Type = BOOL_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				BoolType, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  blitz::Array<unsigned char,1> wval(v.extent(0));
  for (int i=0; i<v.extent(0); i++)
    wval(i) = v(i) ? (unsigned char)1 : (unsigned char)0;
  status = H5Dwrite(newVar->DataSetID, BoolType, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, wval.data());
  if (status < 0)
    cerr << "Error writing bool to HDF5 file in WriteVar.\n";


  // Free DSprops
  H5Pclose (DSprops);
}


void IOTreeHDF5Class::WriteVar(string name, blitz::Array<bool,2> &v)
{
  const int RANK=2;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->BoolType = BoolType;
  newVar->Name = name;
  newVar->Type = BOOL_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				BoolType, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  blitz::Array<unsigned char,2> wval(v.extent(0),v.extent(1));
  for (int i=0; i<v.extent(0); i++)
    for (int j=0; j<v.extent(1); j++)
      wval(i,j) = v(i,j) ? (unsigned char)1 : (unsigned char)0;
  status = H5Dwrite(newVar->DataSetID, BoolType, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, wval.data());
  if (status < 0)
    cerr << "Error writing bool to HDF5 file in WriteVar.\n";

  // Free DSprops
  H5Pclose (DSprops);
}



void IOTreeHDF5Class::WriteVar(string name, blitz::Array<bool,3> &v)
{
  const int RANK=3;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->BoolType = BoolType;
  newVar->Name = name;
  newVar->Type = BOOL_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				BoolType, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  blitz::Array<unsigned char,3> wval(v.extent(0), v.extent(1), v.extent(2));
  for (int i=0; i<v.extent(0); i++)
    for (int j=0; j<v.extent(1); j++)
      for (int k=0; k<v.extent(2); k++)
	wval(i,j,k) = v(i,j,k) ? (unsigned char)1 : (unsigned char)0;
  status = H5Dwrite(newVar->DataSetID, BoolType, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, wval.data());
  if (status < 0)
    cerr << "Error writing bool to HDF5 file in WriteVar.\n";

  // Free DSprops
  H5Pclose (DSprops);
}



void IOTreeHDF5Class::WriteVar(string name, blitz::Array<bool,4> &v)
{
  const int RANK=4;
  /// Create new variable
  VarHDF5Class *newVar = new VarHDF5Class;
  newVar->BoolType = BoolType;
  newVar->Name = name;
  newVar->Type = BOOL_TYPE;
  newVar->Dim = RANK;
  newVar->Dimensions.resize(RANK);
  for (int i=0; i<RANK; i++)
    newVar->Dimensions(i) = v.extent(i);

  /// Add it to the list
  VarList.push_back(newVar);

  // Now do HDF5 stuff
  hid_t dataspace_id, status, DSprops;

  hsize_t dim[RANK], maxdim[RANK];
  for(int i=0; i<RANK; i++)
    dim[i] = v.extent(i);
  // First dimension can be extended indefinitely
  maxdim[0] = H5S_UNLIMITED;
  
  // Rest of dimensions can't be extended
  for(int i=1; i<RANK; i++)
    maxdim[i] = v.extent(i);

  newVar->DataSpaceID = H5Screate_simple(RANK, dim, maxdim);
  DSprops = H5Pcreate (H5P_DATASET_CREATE);
  /// Chunk_dims tells us how big of chunks to allocate in the file at
  /// a time.
  hsize_t chunk_dims[RANK];
  chunk_dims[0]=30;
  for (int i=1;i<RANK;i++)
    chunk_dims[i]=v.extent(i);

  status = H5Pset_chunk(DSprops, RANK, chunk_dims);  
  assert(status>=0);
 
  // Actually create the dataspace
  newVar->DataSetID = H5Dcreate(GroupID, name.c_str(),
				BoolType, newVar->DataSpaceID,
				DSprops);
  ///Extend dataset to size of v
  status=H5Dextend(newVar->DataSetID,dim); 
  assert(status>=0);
  // Do the writing.
  blitz::Array<unsigned char,4> wval(v.extent(0), v.extent(1), 
			      v.extent(2), v.extent(3));
  for (int i=0; i<v.extent(0); i++)
    for (int j=0; j<v.extent(1); j++)
      for (int k=0; k<v.extent(2); k++)
	for (int l=0; l<v.extent(3); l++)
	  wval(i,j,k,l) = v(i,j,k,l) ? (unsigned char)1 : (unsigned char)0;
  status = H5Dwrite(newVar->DataSetID, BoolType, 
		    H5S_ALL, H5S_ALL, H5P_DEFAULT, wval.data());
  if (status < 0)
    cerr << "Error writing int to HDF5 file in WriteVar.\n";

  // Free DSprops
  H5Pclose (DSprops);
}



