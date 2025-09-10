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
  void Atom_writer::dump_xyz_mpi_shared_atoms(int64_t, double)
  {
#if defined(CAVIAR_WITH_MPI)

    auto &pos = atom_data->atom_struct_owned.position;
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &acc = atom_data->atom_struct_owned.acceleration;
    const auto &id = atom_data->atom_struct_owned.id;
    //const auto &type = atom_data->atom_struct_owned.type;

    unsigned int pos_size = pos.size(); // num_local_atoms;
    const unsigned nprocs = comm->nprocs;
    //unsigned int num_total_atoms = pos.size();

    unsigned int num_active_atoms = 0;
    for (unsigned int i = 0; i < pos_size; ++i)
    {    
      if ( atom_data->atom_struct_owned.mpi_rank[i] == my_mpi_rank) num_active_atoms++;
    }
    unsigned int send_num = num_active_atoms;

    std::vector<unsigned int> recv_num(nprocs, 0);
    //==================================// num_local_atoms send
    if (my_mpi_rank != 0)
    {

      MPI_Send(&send_num, 1, MPI::UNSIGNED, 0, 0, MPI_COMM_WORLD);
    }
    if (my_mpi_rank == 0)
    {
      for (unsigned i = 1; i < nprocs; ++i)
      {
        MPI_Recv(&recv_num[i], 1, MPI::UNSIGNED, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //num_total_atoms += recv_num[i];
      }
    }
    auto &mpinf = mpi_packet_info_owned;

    if (!mpinf.initialized)
    {
      mpinf.initialized = true;
      mpinf.total = 0;
      //if (output_id)
      {
        mpinf.id = mpinf.total;
        mpinf.total += 1;
      }
      //mpinf.type = mpinf.total;
      //mpinf.total += 1;

      mpinf.pos = mpinf.total;
      mpinf.total += 3;

      if (xyz_output_velocity)
      {
        mpinf.vel = mpinf.total;
        mpinf.total += 3;
      }
      if (xyz_output_acceleration)
      {
        mpinf.acc = mpinf.total;
        mpinf.total += 3;
      }
    }

    std::vector<double> send_data;              // rank != 0
    std::vector<std::vector<double>> recv_data; // rank 0

    if (my_mpi_rank != 0)
    {

      send_data.resize(mpinf.total * send_num, 0);
    }

    if (my_mpi_rank == 0)
    {
      recv_data.resize(nprocs);
      for (unsigned i = 1; i < nprocs; ++i)
      {
        if (recv_num[i] > 0)
          recv_data[i].resize(mpinf.total * recv_num[i], 0);
      }
    }

    if (my_mpi_rank != 0)
    {
      auto N = send_num;

      int i = 0;
      for (unsigned int j = 0; j < pos_size; ++j)
      {
        if ( atom_data->atom_struct_owned.mpi_rank[j] != my_mpi_rank) continue;

        //if (output_id)
          send_data[mpinf.id * N + i] = id[j];

        //send_data[mpinf.type * N + i] = type[j];

        send_data[(mpinf.pos * N) + (3 * i) + 0] = pos[j].x;
        send_data[(mpinf.pos * N) + (3 * i) + 1] = pos[j].y;
        send_data[(mpinf.pos * N) + (3 * i) + 2] = pos[j].z;

        if (xyz_output_velocity)
        {
          send_data[(mpinf.vel * N) + (3 * i) + 0] = vel[j].x;
          send_data[(mpinf.vel * N) + (3 * i) + 1] = vel[j].y;
          send_data[(mpinf.vel * N) + (3 * i) + 2] = vel[j].z;
        }

        if (xyz_output_acceleration)
        {
          send_data[(mpinf.acc * N) + (3 * i) + 0] = acc[j].x;
          send_data[(mpinf.acc * N) + (3 * i) + 1] = acc[j].y;
          send_data[(mpinf.acc * N) + (3 * i) + 2] = acc[j].z;
        }
        ++i;
      }
    }
    //-----------------------------------------------//
    if (my_mpi_rank != 0)
      if (send_num > 0)

        MPI_Send(send_data.data(), mpinf.total * send_num, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // TAG 1

    if (my_mpi_rank == 0)
      for (unsigned i = 1; i < nprocs; ++i)
        if (recv_num[i] > 0)

          MPI_Recv(recv_data[i].data(), mpinf.total * recv_num[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // TAG 1

    //-----------------------------------------------//
    if (my_mpi_rank != 0)
      return;

    Vector3d<double> p_o{0, 0, 0};

    if (position_offset != nullptr)
      p_o = position_offset->current_value;

    
    // sync all atom data in mpi_rank == 0
    for (unsigned i = 1; i < nprocs; ++i)
    {
      int N = recv_num[i];

      for (unsigned int j = 0; j < recv_num[i]; ++j)
      {

        //auto type_j = recv_data[i][mpinf.type * N + j];
        auto id_j = recv_data[i][mpinf.id * N + j];
        auto k = atom_data->atom_id_to_index[id_j];
        pos[k].x = recv_data[i][(mpinf.pos * N) + (3 * j) + 0];
        pos[k].y = recv_data[i][(mpinf.pos * N) + (3 * j) + 1];
        pos[k].z = recv_data[i][(mpinf.pos * N) + (3 * j) + 2];

        if (xyz_output_velocity)
        {
          vel[k].x = recv_data[i][(mpinf.vel * N) + (3 * j) + 0];
          vel[k].y = recv_data[i][(mpinf.vel * N) + (3 * j) + 1];
          vel[k].z = recv_data[i][(mpinf.vel * N) + (3 * j) + 2];
        }

        if (xyz_output_acceleration)
        {
          acc[k].x = recv_data[i][(mpinf.acc * N) + (3 * j) + 0];
          acc[k].y = recv_data[i][(mpinf.acc * N) + (3 * j) + 1];
          acc[k].z = recv_data[i][(mpinf.acc * N) + (3 * j) + 2];
        }
      }
    }
    dump_xyz_serial(0,false);

/*
    ofs_xyz << num_total_atoms << "\nAtom\n";

    for (unsigned int i = 0; i < pos_size; ++i)
    {

      if ( atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank) continue;

      ofs_xyz << type[i];
      if (output_id)
      {

        ofs_xyz << " " << id[i];
      }

      ofs_xyz << " " << pos[i].x + p_o.x << " " << pos[i].y + p_o.y << " " << pos[i].z + p_o.z;

      if (xyz_output_velocity)
        ofs_xyz << " " << vel[i].x << " " << vel[i].y << " " << vel[i].z;
      if (xyz_output_acceleration)
        ofs_xyz << " " << acc[i].x << " " << acc[i].y << " " << acc[i].z;
      ofs_xyz << "\n";
    }

    for (unsigned i = 1; i < nprocs; ++i)
    {

      int N = recv_num[i];

      for (unsigned int j = 0; j < recv_num[i]; ++j)
      {

        auto type_j = recv_data[i][mpinf.type * N + j];

        auto pos_x = recv_data[i][(mpinf.pos * N) + (3 * j) + 0];
        auto pos_y = recv_data[i][(mpinf.pos * N) + (3 * j) + 1];
        auto pos_z = recv_data[i][(mpinf.pos * N) + (3 * j) + 2];

        ofs_xyz << type_j;

        if (output_id)
        {
          auto id_j = recv_data[i][mpinf.id * N + j];

          ofs_xyz << " " << id_j;
        }

        ofs_xyz << " " << pos_x + p_o.x << " " << pos_y + p_o.y << " " << pos_z + p_o.z;

        if (xyz_output_velocity)
        {
          auto vel_x = recv_data[i][(mpinf.vel * N) + (3 * j) + 0];
          auto vel_y = recv_data[i][(mpinf.vel * N) + (3 * j) + 1];
          auto vel_z = recv_data[i][(mpinf.vel * N) + (3 * j) + 2];
          ofs_xyz << " " << vel_x << " " << vel_y << " " << vel_z;
        }

        if (xyz_output_acceleration)
        {
          auto acc_x = recv_data[i][(mpinf.acc * N) + (3 * j) + 0];
          auto acc_y = recv_data[i][(mpinf.acc * N) + (3 * j) + 1];
          auto acc_z = recv_data[i][(mpinf.acc * N) + (3 * j) + 2];
          ofs_xyz << " " << acc_x << " " << acc_y << " " << acc_z;
        }
        ofs_xyz << "\n";
      }
    }

    ofs_xyz << std::flush;
*/
#endif
  }

    //================================================
  //                                              ||
  //================================================
  void Atom_writer::dump_xyz_ghost_mpi_shared_atoms(int64_t, double)
  {
    //dump_xyz_ghost_mpi_shared_atoms(int64_t, double);
  }

} // writer

}
