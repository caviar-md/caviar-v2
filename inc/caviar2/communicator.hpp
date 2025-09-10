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

#pragma once

#include <string>

#if defined(CAVIAR_WITH_MPI)
// MPI_Datatype mpi_fc_vector_type;
#endif

namespace caviar2
{
  class Caviar2;

  /**
   * This class handles MPI process base communications.
   *
   */
  class Communicator
  {
  public:
    Communicator(Caviar2 *);
    Communicator();
    virtual void set_caviar(Caviar2 *c);

    // broadcast a variable from the root process to others
    void broadcast(bool &);
    void broadcast(size_t &);
    void broadcast(size_t &, char *);
    void broadcast(std::string &);

#if defined(CAVIAR_WITH_MPI)
    // MPI_Datatype mpi_fc_vector_type;
#endif

    int me, nprocs; // MPI process rank and number of processes
    Caviar2 *caviar_ = nullptr;
  };
}
