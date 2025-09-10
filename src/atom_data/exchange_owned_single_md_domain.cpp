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

#include <algorithm>

namespace caviar2 {

//======================================================
//                                                    ||
//                                                    ||
//======================================================

bool Atom_data::exchange_owned_single_md_domain(long) // timestep
{
  if (domain == nullptr)
    caviar_->log.error_all("Atom_data::exchange_owned: domain = nullptr");

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
  if (domain->me != 0)
    return false;
#endif

  bool update_verlet_list = false;

  auto &pos = atom_struct_owned.position;

  auto pos_size = pos.size();

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < pos_size; ++i)
  {
    caviar2::Vector3d<int> msd {0,0,0};
    pos[i] = domain->fix_position(pos[i], msd, update_verlet_list);

    if (msd_process)
      atom_struct_owned.msd_domain_cross[i] +=  msd;
  }
  return update_verlet_list;
}
//======================================================
//                                                    ||
//                                                    ||
//======================================================

}
