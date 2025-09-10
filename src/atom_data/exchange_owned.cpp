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

#include "caviar2/atom_data.hpp"
#include "caviar2/communicator.hpp"
#include "caviar2/caviar2.hpp"
#include "caviar2/domain.hpp"


namespace caviar2 {


bool Atom_data::exchange_owned(long i) // timestep
{
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  return exchange_owned_single_md_domain();

#elif defined(CAVIAR_WITH_MPI)

  return exchange_owned_mpi_shared_atoms(i);
  // return exchange_owned_mpi(i);

#else

  return exchange_owned_single_md_domain(i);

#endif
}

}
