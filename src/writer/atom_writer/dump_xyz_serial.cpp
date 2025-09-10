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

#include "caviar2/writer/atom_writer.hpp"
#include "caviar2/atom_data.hpp"

#include "caviar2/communicator.hpp"
#include "caviar2/time_utility.hpp"
#include "caviar2/tools/time_function_3d.hpp"
#ifdef CAVIAR_WITH_MPI
#include <mpi.h>
#endif
namespace caviar2 {

namespace writer
{

  //================================================
  //                                              ||
  //================================================
  void Atom_writer::dump_xyz_serial(int64_t, double)
  {
    
 
    auto &pos = atom_data->atom_struct_owned.position;
    auto &type = atom_data->atom_struct_owned.type;
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &acc = atom_data->atom_struct_owned.acceleration;
    auto &id = atom_data->atom_struct_owned.id;

    auto nta = atom_data->atom_struct_owned.position.size();

    Vector3d<double> p_o{0, 0, 0};
    if (position_offset != nullptr)
      p_o = position_offset->current_value;
#if defined(CAVIAR_WITH_MPI)
    if (xyz_mpi_rank0)
#endif
    {
      ofs_xyz << nta << "\nAtom\n";

      for (unsigned int i = 0; i < nta; ++i)
      {
        ofs_xyz << type[i];
        if (xyz_output_id)
          ofs_xyz << " " << id[i];
        ofs_xyz << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
        if (xyz_output_velocity)
          ofs_xyz << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
        if (xyz_output_acceleration)
          ofs_xyz << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
        ofs_xyz << "\n";
      }

      ofs_xyz << std::flush;
    }
  }

  void Atom_writer::dump_xyz_mpi_per_process(int64_t, double)
  {

#if defined(CAVIAR_WITH_MPI)

    auto &pos = atom_data->atom_struct_owned.position;
    auto &type = atom_data->atom_struct_owned.type;
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &acc = atom_data->atom_struct_owned.acceleration;
    auto &id = atom_data->atom_struct_owned.id;

    auto nta = atom_data->atom_struct_owned.position.size();

    Vector3d<double> p_o{0, 0, 0};
    if (position_offset != nullptr)
      p_o = position_offset->current_value;
    if (xyz_mpi_per_process)
    {
      unsigned int num_active_atoms = 0;
      for (unsigned int i = 0; i < nta; ++i)
      {
        if (atom_data->atom_struct_owned.mpi_rank[i] == my_mpi_rank)
          num_active_atoms++;
      }

      ofs_xyz_mpi << num_active_atoms << "\nAtom\n";

      for (unsigned int i = 0; i < nta; ++i)
      {
        if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
          continue;

        ofs_xyz_mpi << type[i];
        if (xyz_output_id)
          ofs_xyz_mpi << " " << id[i];

        ofs_xyz_mpi << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
        if (xyz_output_velocity)
          ofs_xyz_mpi << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
        if (xyz_output_acceleration)
          ofs_xyz_mpi << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
        ofs_xyz_mpi << "\n";
      }

      ofs_xyz_mpi << std::flush;
    }
#endif

  }

  //================================================
  //                                              ||
  //================================================
  void Atom_writer::dump_xyz_ghost_serial(int64_t, double)
  {
    auto &pos = atom_data->atom_struct_ghost.position;
    auto &type = atom_data->atom_struct_ghost.type;
    auto &vel = atom_data->atom_struct_ghost.velocity;
    auto &acc = atom_data->atom_struct_ghost.acceleration;
    auto &id = atom_data->atom_struct_ghost.id;
    auto nta = atom_data->atom_struct_ghost.position.size();

    Vector3d<double> p_o{0, 0, 0};
    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    if (xyz_ghost_mpi_rank0)
    {
      ofs_xyz_ghost << nta << "\nAtom\n";

      for (unsigned int i = 0; i < nta; ++i)
      {
        ofs_xyz_ghost << type[i];
        if (xyz_output_id)
          ofs_xyz_ghost << " " << id[i];
        ofs_xyz_ghost << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
        if (xyz_output_velocity)
          ofs_xyz_ghost << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
        if (xyz_output_acceleration)
          ofs_xyz_ghost << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
        ofs_xyz_ghost << "\n";
      }

      ofs_xyz_ghost << std::flush;
    }
  }
  void Atom_writer::dump_xyz_ghost_mpi_per_process(int64_t, double)
  {
    auto &pos = atom_data->atom_struct_ghost.position;
    auto &type = atom_data->atom_struct_ghost.type;
    auto &vel = atom_data->atom_struct_ghost.velocity;
    auto &acc = atom_data->atom_struct_ghost.acceleration;
    auto &id = atom_data->atom_struct_ghost.id;
    auto nta = atom_data->atom_struct_ghost.position.size();

    Vector3d<double> p_o{0, 0, 0};
    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    if (xyz_ghost_mpi_per_process)
    {
      ofs_xyz_ghost_mpi << nta << "\nAtom\n";

      for (unsigned int i = 0; i < nta; ++i)
      {
        ofs_xyz_ghost_mpi << type[i];
        if (xyz_output_id)
          ofs_xyz_ghost_mpi << " " << id[i];

        ofs_xyz_ghost_mpi << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;
        if (xyz_output_velocity)
          ofs_xyz_ghost_mpi << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
        if (xyz_output_acceleration)
          ofs_xyz_ghost_mpi << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
        ofs_xyz_ghost_mpi << "\n";
      }

      ofs_xyz_ghost_mpi << std::flush;
    }
  }

} // writer

}
