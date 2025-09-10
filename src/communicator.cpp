/*
 * -----------------------------------------------------------------------------
 * CAVIAR2 - C++ Library
 *
 * Copyright (c) 2025 Morad Biagooi and Ehsan Nedaaee Oskoee
 * All rights reserved.
 *
 * License: To be determined.
 * This file is provided "as is", without warranty of any kind.
 * You may not distribute this code until a license is finalized.
 *
 * -----------------------------------------------------------------------------
 */

#include "caviar2/communicator.hpp"
#include "caviar2/vector3d.hpp"
#if defined(CAVIAR_WITH_MPI)
#include <mpi.h>
#endif

namespace caviar2
{
  Communicator::Communicator() {}

  Communicator::Communicator(Caviar2 *fptr) : caviar_{fptr}
  {
#if defined(CAVIAR_WITH_MPI)

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /*
    // Adding a MPI type for caviar2::Vector
    const int nitems = 3;
    int blocklengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    // MPI_Datatype mpi_fc_vector_type;
    MPI_Aint offsets[3];

    offsets[0] = offsetof(caviar2::Vector3d<double>, x);
    offsets[1] = offsetof(caviar2::Vector3d<double>, y);
    offsets[2] = offsetof(caviar2::Vector3d<double>, z);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_fc_vector_type);
    MPI_Type_commit(&mpi_fc_vector_type);
    */
#endif
  }

  void Communicator::set_caviar(caviar2::Caviar2 *c)
  {
    caviar_ = c;
  }

  void Communicator::broadcast(bool &flag)
  {
#if defined(CAVIAR_WITH_MPI)
    MPI_Bcast(&flag, 1, MPI::BOOL, 0, MPI_COMM_WORLD);
#else
    std::cout << "Communicator::broadcast " << flag << std::endl;
#endif
  }

  void Communicator::broadcast(size_t &n)
  {
#if defined(CAVIAR_WITH_MPI)
    MPI_Bcast(&n, 1, MPI::INT, 0, MPI_COMM_WORLD);
#else
    std::cout << "Communicator::broadcast " << n << std::endl;
#endif
  }

  void Communicator::broadcast(size_t &n, char *str)
  {
#if defined(CAVIAR_WITH_MPI)
    MPI_Bcast(&n, 1, MPI::INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(str, n, MPI::CHAR, 0, MPI_COMM_WORLD);
#else
    std::cout << "Communicator::broadcast " << n << " " << str << std::endl;
#endif
  }

  void Communicator::broadcast(std::string &str)
  {
#if defined(CAVIAR_WITH_MPI)

    int n = me == 0 ? str.length() : 0;
    MPI_Bcast(&n, 1, MPI::INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // The reason behind '[n+1]' for 'char *tmp' is because of 'strcpy':
    // (from http://www.cplusplus.com/reference/cstring/strcpy/)
    // To avoid overflows, the size of the array pointed by destination shall be
    // long enough to contain the same C string as source (including the
    // terminating null character), and should not overlap in memory with source.
    char *tmp = new char[n + 1];

    strcpy(tmp, str.c_str());

    MPI_Bcast(tmp, n, MPI::CHAR, 0, MPI_COMM_WORLD);

    if (tmp)
    {
      if (tmp[0])
      {
        str.assign(const_cast<const char *>(tmp), n);
      }
      else
      {
      }
    }
    else
    {
    }

    MPI_Barrier(MPI_COMM_WORLD);

    delete[] tmp;

#else
    std::cout << "Communicator::broadcast is called in non-mpi mode" << str << std::endl;
#endif
  }

}
