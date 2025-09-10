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

void Atom_data::synch_owned_data(int)
{
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

  auto mpi_fc_vector_type = comm->mpi_fc_vector_type;
  auto my_mpi_rank = comm->me;
  if (comm->nprocs == 1)
    return;

  int root_rank = 0;

  int pos_size = atom_struct_owned.position.size();
  int root_pos_size = pos_size;

  int type_size = atom_struct_owned.type.size();
  // int root_type_size = type_size;

  // int mass_size = atom_type_params.mass.size();
  // int root_mass_size = mass_size;

  // int charge_size = atom_type_params.charge.size();
  // int root_charge_size = charge_size;

  MPI_Bcast(&root_pos_size, 1, MPI::INT, 0, mpi_comm);

  // XXX not necessary yet.
  // MPI_Bcast (&synch_owned_data_bcast_details,   1, MPI::BOOL, 0, mpi_comm);

  // if (synch_owned_data_bcast_details)
  // {
  //   MPI_Bcast(&root_type_size, 1, MPI::INT, 0, mpi_comm);
  //   MPI_Bcast(&root_mass_size, 1, MPI::INT, 0, mpi_comm);
  //   MPI_Bcast(&root_charge_size, 1, MPI::INT, 0, mpi_comm);
  // }

  if (my_mpi_rank != root_rank)
  {

    if (pos_size != root_pos_size)
    {
      atom_struct_owned.position.resize(root_pos_size, caviar2::Vector3d<double>{0, 0, 0});
      atom_struct_owned.velocity.resize(root_pos_size, caviar2::Vector3d<double>{0, 0, 0});
      atom_struct_owned.acceleration.resize(root_pos_size, caviar2::Vector3d<double>{0, 0, 0});
    }

    if (type_size != root_pos_size)
    {
      atom_struct_owned.type.resize(root_type_size, 0);
    }

    // if (mass_size != root_mass_size)
    // {
    //   atom_type_params.mass.resize(root_mass_size, 1.0);
    // }

    // if (charge_size != root_charge_size)
    // {
    //   atom_type_params.charge.resize(root_charge_size, 0.0);
    // }
  }

  MPI_Bcast(&atom_struct_owned.position[0], root_pos_size, mpi_fc_vector_type, 0, mpi_comm);

  // if (synch_owned_data_bcast_details)
  // {
  //   MPI_Bcast(&atom_struct_owned.type[0], root_type_size, MPI::INT, 0, mpi_comm);

  //   MPI_Bcast(&atom_type_params.mass[0], root_mass_size, MPI::DOUBLE, 0, mpi_comm);

  //   MPI_Bcast(&atom_type_params.charge[0], root_charge_size, MPI::DOUBLE, 0, mpi_comm);
  // }

  synch_owned_data_bcast_details = false;

#endif
}

}
