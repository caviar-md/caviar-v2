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
#include "caviar2/domain.hpp"


namespace caviar2 {


void Atom_data::exchange_ghost(long i) // timestep
{
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  exchange_ghost_single_md_domain();

#elif defined(CAVIAR_WITH_MPI)

  exchange_ghost_mpi_shared_atoms(i);
  // exchange_ghost_mpi(i);

#else

  exchange_ghost_single_md_domain(i);

#endif
}

//  if (self_ghost_check())
//    caviar_->log.error_all (FC_FILE_LINE_FUNC_PARSE, "Self ghost can happen. Force field cutoff is larger than half of a domain.");
/*
bool ::self_ghost_check () {
  const auto x_llow = domain->lower_local.x;
  const auto x_lupp = domain->upper_local.x;
  const auto y_llow = domain->lower_local.y;
  const auto y_lupp = domain->upper_local.y;
  const auto z_llow = domain->lower_local.z;
  const auto z_lupp = domain->upper_local.z;

  const auto x_width = x_lupp - x_llow;
  const auto y_width = y_lupp - y_llow;
  const auto z_width = z_lupp - z_llow;

  const auto cutoff = force_field->cutoff;
  if (2*cutoff>x_width || 2*cutoff>y_width || 2*cutoff>z_width)
    return true;
  return false;
}*/

}
