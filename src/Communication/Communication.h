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

#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "../Blitz.h"
#include <fstream>

extern ostream perr;

namespace COMM
{
  using blitz::Array;
#ifdef USE_MPI
  inline void Init (int argc, char **argv)
  {
    MPI_Init(&argc, &argv);
    int proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    if (proc == 0)
      perr.rdbuf(cerr.rdbuf());
//      else 
//        perr.open("/dev/null", ios::out); 
  }
  inline void Finalize ()
  {
    MPI_Finalize();
  }
  inline int WorldProc()
  {
    int proc;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    return (proc);
  }
#else // Serial version
  inline void Init (int argc, char **argv)
  {
    perr.rdbuf(cerr.rdbuf());
  }
  inline void Finalize ()
  {
  }
  inline int WorldProc()
  {
    return (0);
  }

#endif
}  


class CommStatusClass
{
public:
  int Source;
  int Tag;
  int Error;
  int Length;
};


class CommunicatorClass
{
public:
#ifdef USE_MPI
  ///If we are in parallel mode, we need an MPI communicator
  MPI_Comm MPIComm;  
  
  /// Sets this communicator to be that of all the processes
  /// (i.e. MPI_WORLD)
  void SetWorld();
  int MyProc();
  int NumProcs();
  string MyHost();
  void Send (void *sendBuf, int count, MPI_Datatype datatype,
	     int dest, int tag);
  void Send (int toProc, blitz::Array<double,1> &buff);
  void Send (int toProc, blitz::Array<int,1> &buff);

  void Send    (int   toProc, int val);
  void Receive (int fromProc, int &val);
  void Send    (int   toProc, double val);
  void Receive (int fromProc, double &val);

  void Broadcast (int root, int &val);
  void Broadcast (int root, bool &val);
  void Broadcast (int root, double &val);
  void Broadcast (int root, blitz::Array<int,1> &buff);
  void Broadcast (int root, blitz::Array<double,1> &buff);
  void Broadcast (int root, blitz::Array<double,2> &buff);
  void Broadcast (int root, blitz::Array<double,3> &buff);
  void Broadcast (int root, blitz::Array<complex<double>,1> &buff);
  void Broadcast (int root, blitz::Array<complex<double>,2> &buff);
  void Broadcast (int root, blitz::Array<Vec1,1> &buff);
  void Broadcast (int root, blitz::Array<Vec2,1> &buff);
  void Broadcast (int root, blitz::Array<Vec3,1> &buff);
  void Broadcast (int root, blitz::Array<double,4> &buff);
  void Broadcast (int root, blitz::Array<complex<double>,4> &buff);
  void Receive (void *recvBuf, int count, MPI_Datatype datatype,
		int source, int tag);
  void Receive (int toProc, blitz::Array<double,1> &buff);
  void Receive (int toProc, blitz::Array<int,1> &buff);
  bool Probe(int source, int tag, CommStatusClass &status);
  void GatherProd(double &send, double &recv, int root);
  void Gather (blitz::Array<complex<double>,1> &sendVec, 
	       blitz::Array<complex<double>,1> &recvVec, 
	       blitz::Array<int,1> &recvCounts, int root=0);
  void Gather (blitz::Array<TinyVector<double,3>,1> &sendVec, 
	       blitz::Array<TinyVector<double,3>,1> &recvVec, 
	       blitz::Array<int,1> &recvCounts, int root=0);
  void Gather (blitz::Array<double,1> &sendVec, blitz::Array<double,2> &recvMat,
	       int root=0);
  void Gather (blitz::Array<int,1> &sendVec, blitz::Array<int,2> &recvMat,
	       int root=0);
  void AllGather(void *sendbuf, int sendcount, 
		 MPI_Datatype sendtype, 
		       void* recvbuf, int recvcount,
		 MPI_Datatype recvtype);
  void AllGather (blitz::Array<double,1> &SendVec, 
		  blitz::Array<double,1> &RecvVec);
  void AllGather (blitz::Array<int,1> &SendVec, 
		  blitz::Array<int,1> &RecvVec);
  /// This function gathers all rows of a matrix to all processors.
  /// It assumes that each processor receives an equal number of rows,
  /// with any left over being distributed to the low number
  /// processors. E.g. If there are 8 rows and 3 processors, procs 0
  /// and 1 would get 3 rows and proc 2 would get 2 rows.
  void AllGatherRows (blitz::Array<complex<double>,2> &mat);
  void AllGatherRows (blitz::Array<double,2> &mat);

  /// This function uses the same division stragegy as above, but
  /// gathers single elements, instead.
  void AllGatherVec (blitz::Array<double,1> &vec);
  void AllGatherVec(blitz::Array<int,1> &vec);


  void Split (int color, CommunicatorClass &newComm);
  void Subset (blitz::Array<int,1> &ranks, CommunicatorClass &newComm);

//   ///Sends and receives an array of dVec
//   void SendReceive (int sendProc, const blitz::Array<Vec3,1> &sendBuff,
// 		    int recvProc,       blitz::Array<Vec3,1> &recvBuff);

//   ///Sends and receives an array of dVec
//   void SendReceive (int sendProc, const blitz::Array<Vec2,1> &sendBuff,
//		    int recvProc,       blitz::Array<Vec2,1> &recvBuff);
  


  template<int N>
  void Send (int toProc, TinyVector<double,N> val)
  {
    MPI_Send(&val[0], N, MPI_DOUBLE, toProc, 1, MPIComm);
  }

  template<int N>
  void Receive (int fromProc, TinyVector<double,N> &val)
  {
    MPI_Status status;
    MPI_Recv(&val[0], N, MPI_DOUBLE, fromProc, 1, MPIComm, &status);
  }


  template<int N>
  void Sum (Array<double,N> &sendBuff, Array<double,N> &recvBuff)
  {
    double *sendPtr = sendBuff.data();
    double *recvPtr = recvBuff.data();
    int count = sendBuff.size();
    MPI_Reduce(sendPtr, recvPtr, count, MPI_DOUBLE, MPI_SUM, 0, 
	       MPIComm);
  }




  template<int N>
  void Sum (blitz::Array<complex<double>,N> Ain, blitz::Array<complex<double>,N> Aout)
  {
    assert (Ain.size() == Aout.size());
    MPI_Reduce ((void*)Ain.data(), (void*)Aout.data(), 
		2*Ain.size(), MPI_DOUBLE, MPI_SUM, 0, MPIComm);
  }

  template<int N>
  void AllSum (blitz::Array<complex<double>,N> Ain, blitz::Array<complex<double>,N> Aout)
  {
    assert (Ain.size() == Aout.size());
    MPI_Allreduce ((void*)Ain.data(), (void*)Aout.data(), 
		   2*Ain.size(), MPI_DOUBLE, MPI_SUM, MPIComm);
  }

  template<int N>
  void Broadcast (int root, TinyVector<double,N> &vec)
  {
    MPI_Bcast(&(vec[0]), N, MPI_DOUBLE, root, MPIComm);
  }

  template<int N, int M>
  void Broadcast (int root, TinyMatrix<double,N,M> &mat)
  {
    MPI_Bcast(&(mat(0,0)), N*M, MPI_DOUBLE, root, MPIComm);
  }


  template<int N>
  void Broadcast (int root, Array<int,N> &a)
  {
    assert (a.isStorageContiguous());
    MPI_Bcast(a.data(), a.size(), MPI_INT, root, MPIComm);
  }

  template<int N>
  void Broadcast (int root, Array<TinyVector<int,N>,1> &vec)
  {
    MPI_Bcast(&(vec(0)[0]), N*vec.size(), MPI_INT, root, MPIComm);
  }


  template<int N>
  void SendReceive (int sendProc, const Array<double,N> &sendBuff,
		    int recvProc,       Array<double,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), sendBuff.size(), MPI_DOUBLE, sendProc, 2,
		 (void*)recvBuff.data(), recvBuff.size(), MPI_DOUBLE, recvProc, 2,
		 MPIComm, &status);
  }

  template<int N, int M>
  void SendReceive (int sendProc, const blitz::Array<TinyVector<double,M>,N> &sendBuff,
		    int recvProc,       blitz::Array<TinyVector<double,M>,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), M*sendBuff.size(), MPI_DOUBLE, sendProc, 2,
		 (void*)recvBuff.data(), M*recvBuff.size(), MPI_DOUBLE, recvProc, 2,
		 MPIComm, &status);
  }
  
  template<int N>
  void SendReceive (int sendProc, const blitz::Array<complex<double>,N> sendBuff,
		    int recvProc,       blitz::Array<complex<double>,N> recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), 2*sendBuff.size(), MPI_DOUBLE, sendProc, 2,
		 (void*)recvBuff.data(), 2*recvBuff.size(), MPI_DOUBLE, recvProc, 2,
		 MPIComm, &status);
  }

  template<int N>
  void Send (int sendProc, const blitz::Array<complex<double>,N> sendBuff)
  {
    MPI_Send((void*)sendBuff.data(), 2*sendBuff.size(), MPI_DOUBLE, sendProc, 2,
	     MPIComm);
  }
  
  template<int N>
  void Receive (int recvProc, blitz::Array<complex<double>,N> recvBuff)
  {
    MPI_Status status;
    MPI_Recv((void*)recvBuff.data(), 2*recvBuff.size(), MPI_DOUBLE, recvProc, 2,
	     MPIComm, &status);
  }

  template<int N, int M>
  void SendReceive (int sendProc, const blitz::Array<TinyVector<complex<double>,M>,N> &sendBuff,
		    int recvProc,       blitz::Array<TinyVector<complex<double>,M>,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), 2*M*sendBuff.size(), MPI_DOUBLE, sendProc, 2,
		 (void*)recvBuff.data(), 2*M*recvBuff.size(), MPI_DOUBLE, recvProc, 2,
		 MPIComm, &status);
  }

  template<int N>
  void SendReceive (int sendProc, const blitz::Array<int,N> &sendBuff,
		    int recvProc,       blitz::Array<int,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), sendBuff.size(), MPI_INT, sendProc, 3,
		 (void*)recvBuff.data(), recvBuff.size(), MPI_INT, recvProc, 3,
		 MPIComm, &status);
  }

  template<int N, int M>
  void SendReceive (int sendProc, const blitz::Array<TinyVector<int,M>,N> &sendBuff,
		    int recvProc,       blitz::Array<TinyVector<int,M>,N> &recvBuff)
  {
    MPI_Status status;
    MPI_Sendrecv((void*)sendBuff.data(), M*sendBuff.size(), MPI_INT, sendProc, 3,
		 (void*)recvBuff.data(), M*recvBuff.size(), MPI_INT, recvProc, 3,
		 MPIComm, &status);
  }

//   ///Sends and receives an array of double
//   void SendReceive (int sendProc, const blitz::Array<double,1> &sendBuff,
// 		    int recvProc,       blitz::Array<double,1> &recvBuff);

//   ///Sends and receives an array of double
//   void SendReceive (int sendProc, const blitz::Array<double,2> &sendBuff,
// 		    int recvProc,       blitz::Array<double,2> &recvBuff);

//   ///Sends and receives an array of double
//   void SendReceive (int sendProc, const blitz::Array<double,3> &sendBuff,
// 		    int recvProc,       blitz::Array<double,3> &recvBuff);


//   ///Sends and receives an array of complex
//   void SendReceive (int sendProc, const blitz::Array<complex<double>,1> &sendBuff,
// 		    int recvProc,       blitz::Array<complex<double>,1> &recvBuff);
  
//   ///Sends and receives an array of int
//   void SendReceive (int sendProc, const blitz::Array<int,1> &sendBuff,
// 		    int recvProc,       blitz::Array<int,1> &recvBuff);



  ///Sums up the vectors in sendBuff.  Processor 0 only gets the
  ///resulting sum.
  void Sum (blitz::Array<int,1> &sendBuff, blitz::Array<int,1> &recvBuff);
  void Sum (blitz::Array<int,2> &sendBuff, blitz::Array<int,2> &recvBuff);

  ///Sums up the vectors in sendBuff.  Processor 0 only gets the
  ///resulting sum.
  //  void Sum (blitz::Array<double,1> &sendBuff, blitz::Array<double,1> &recvBuff);


  void Sum (blitz::Array<Vec2,1> &sendBuff, blitz::Array<Vec2,1> &recvBuff);
  void Sum (blitz::Array<Vec3,1> &sendBuff, blitz::Array<Vec3,1> &recvBuff);
  ///Sums up all values a.  Only processor 0 gets the result.  All
  ///other processors return 0;
  double Sum (double a);

  /// Sums up all values of a on all processors.  All processors
  ///  get result.
  double AllSum (double a);
  void AllSum (blitz::Array<double,1> &in, blitz::Array<double,1> &out);
  void AllSum (blitz::Array<double,2> &in, blitz::Array<double,2> &out);
  void AllSum (blitz::Array<double,3> &in, blitz::Array<double,3> &out);
  void AllSum (blitz::Array<TinyVector<double,2>,1> &in, 
	       blitz::Array<TinyVector<double,2>,1> &out);
  void AllSum (blitz::Array<TinyVector<double,3>,1> &in, 
	       blitz::Array<TinyVector<double,3>,1> &out);
  void AllAnd (bool &TorF);
  template<int N> 
  inline void AllMax (TinyVector<int,N> &vec) {
    TinyVector<int, N> outVec;
    MPI_Allreduce(&(vec[0]), &(outVec[0]), N, MPI_INT, MPI_MAX, MPIComm);
    vec = outVec;
  }
  void BarrierSync();
  void PrintSync();

  CommunicatorClass()
  {
    SetWorld();
  }
  

#else   // Serial version

  inline void SetWorld()
  {
    // Do nothing
  }
  inline int MyProc()
  {
    return 0;
  }
  inline int NumProcs()
  {
    return 1;
  }

  inline string MyHost()
  {
    return "ThisProc";
  }

  inline void Gather (blitz::Array<complex<double>,1> &sendVec, 
		      blitz::Array<complex<double>,1> &recvVec, 
		      blitz::Array<int,1>& recvCounts, int root=0)
  {
    recvVec = sendVec;
  }
  void GatherProd(double &send, double &recv, int root)
  {
    recv=send;
  }

  inline void Gather (blitz::Array<TinyVector<double,3>,1> &sendVec, 
		      blitz::Array<TinyVector<double,3>,1> &recvVec, 
		      blitz::Array<int,1>& recvCounts, int root=0)
  {
    recvVec = sendVec;
  }
  inline void Gather (blitz::Array<double,1> &sendVec, 
		      blitz::Array<double,2> &recvMat,
		      int root = 0)
  {
    assert (recvMat.rows() == NumProcs());
    assert (recvMat.cols() == sendVec.size());
    recvMat(0,Range::all()) = sendVec;
  }
  inline void Gather (blitz::Array<int,1> &sendVec, 
		      blitz::Array<int,2> &recvMat,
		      int root = 0)
  {
    assert (recvMat.rows() == NumProcs());
    assert (recvMat.cols() == sendVec.size());
    recvMat(0,Range::all()) = sendVec;
  }

  inline void AllGather (blitz::Array<double,1> &SendVec, 
			 blitz::Array<double,1> &RecvVec)
  {
    RecvVec = SendVec;
  }

  inline void AllGather (blitz::Array<int,1> &SendVec, 
			 blitz::Array<int,1> &RecvVec)
  {
    RecvVec = SendVec;
  }

  /// This function gathers all rows of a matrix to all processors.
  /// It assumes that each processor receives an equal number of rows,
  /// with any left over being distributed to the low number
  /// processors. E.g. If there are 8 rows and 3 processors, procs 0
  /// and 1 would get 3 rows and proc 2 would get 2 rows.
  inline void AllGatherRows (blitz::Array<double,2> &mat) 
  {
    // Do nothing
  }
  inline void AllGatherRows (blitz::Array<complex<double>,2> &mat) 
  {
    // Do nothing
  }

  /// This function uses the same division stragegy as above, but
  /// gathers single elements, instead.
  void AllGatherVec (blitz::Array<double,1> &vec) 
  {
    // do nothing
  }

  inline void Split (int color, CommunicatorClass &newComm)
  {
    // Do nothing
  }
  inline void Subset (blitz::Array<int,1> ranks, CommunicatorClass &newComm)
  {
    if (ranks.size() !=1) {
      cerr << "Serial verion of code does not support nontrivial "
	   << "subsets.  Exiting.\n";
      abort();
    }
  }
// <<<<<<< .mine
//   inline void Send (int toProc, blitz::Array<double,1> &buff)
// =======
  template<typename T>
  inline void Send (int toProc, T &val) 
  {
    cerr << "Sends are not support in serial version of CommunicatorClass.\n";
    abort();
  }

  inline void Send (int toProc, blitz::Array<double,1> &buff)
// >>>>>>> .r1224
  {
    cerr << "Sends not supported in serial mode.\n";
    abort();
  }

  inline void Send (int toProc, blitz::Array<int,1> &buff)
  {
    cerr << "Sends not supported in serial mode.\n";
    abort();
  }

  template<int N>
  void SendReceive (int sendProc, const blitz::Array<complex<double>,N> sendBuff,
		    int recvProc,       blitz::Array<complex<double>,N> recvBuff)
  {
    // do nothing for serial version
  }

  template<int N>
  void Send (int sendProc, const blitz::Array<complex<double>,N> sendBuff)
  {
    cerr << "Sends not supported in serial mode.\n";
    abort();
  }
  
  template<typename T>
  void Receive (int recvProc, T &val)
  {
    cerr << "Receives are not supported in serial mode.\n";
    abort();
  }

  template<int N>
  void Receive (int recvProc, blitz::Array<complex<double>,N> recvBuff)
  {
    cerr << "Receives not supported in serial mode.\n";
    abort();
  }

  template<typename T>
  void Broadcast(int root, T &val) { }

  inline void Receive (int toProc, blitz::Array<double,1> &buff)
  { cerr << "Receives not supported in serial mode.\n";  abort(); }
  inline void Receive (int toProc, blitz::Array<int,1> &buff)
  { cerr << "Receives not supported in serial mode.\n";    abort(); }

  template<typename T>
  void SendReceive(int sendProc, const T &sendBuff,
		   int recvProc,       T &recvBuff)
  { recvBuff = sendBuff; }

  ///Sums up all values a.  Only processor 0 gets the result.  All
  ///other processors return 0;
  double Sum (double a)
  { return a; }

  template<class T, int N>
  inline void Sum(blitz::Array<T,N> &sendBuff, blitz::Array<T,N> &recvBuff)
  { recvBuff = sendBuff; }
  
  /// Sums up all values of a on all processors.  All processors get result.
  inline double AllSum (double a)
  {  return a; } 

  template<class T, int N>
  inline void AllSum (blitz::Array<T, N> &in, blitz::Array<T, N> &out)
  { out = in; }

  inline void AllAnd (bool &TorF)
  { /* do nothing */ }

  template<int N> 
  inline void
  AllMax (TinyVector<int,N> &vec) {
    // do nothing
  }


  inline void BarrierSync() {}
  inline void PrintSync() {}

#endif
};  // End class CommunicatorClass
#endif // End ifdef COMMUNICATION_H

