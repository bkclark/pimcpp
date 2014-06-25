class MPI:
    def __init__(self, UseMPI):
        self.UseMPI = UseMPI
        if UseMPI:
            import mpi
            self.mpiLib = mpi
            self.rank = mpi.rank
            self.size = mpi.size
        else:
            self.mpiLib = None
            self.rank = 0
            self.size = 1

    def barrier(self):
        if self.UseMPI:
            self.mpiLib.barrier()

    def allgather(self, a):
        if self.UseMPI:
            return self.mpiLib.allgather(a)
        else:
            return a
