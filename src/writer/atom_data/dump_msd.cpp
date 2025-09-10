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

#include "caviar2/writer/atom_data.hpp"
#include "caviar2/atom_data.hpp"
#include "caviar2/domain.hpp"
#include "caviar2/caviar2.hpp"


namespace caviar2 {

namespace writer
{
  void Atom_data::dump_msd_mpi(int64_t , double )
  {

  }

  void Atom_data::dump_msd_mpi_per_process(int64_t , double )
  {

  }
  
  void Atom_data::dump_msd_mpi_shared_atoms(int64_t , double )
  {

  }

  void Atom_data::dump_msd_serial(int64_t i, double t)
  {

    const auto &pos = atom_data->atom_struct_owned.position;
    const auto &type = atom_data->atom_struct_owned.type;

    auto pos_size = pos.size();
    if (i < msd_initial_step)
      return;
    if (i == msd_initial_step)
    {
      msd_initial_position.resize(pos_size);
      for (unsigned int j = 0; j < pos_size; ++j)
      {
        msd_initial_position[j] = pos[j];
        atom_data->atom_struct_owned.msd_domain_cross[j] = Vector3d<int>{0, 0, 0};
      }
      if (my_mpi_rank == 0)
        ofs_msd << i << " " << t << " "
                << "0.0"
                << "\n";
      return;
    }

    if (pos_size != msd_initial_position.size())
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "  (pos.size != msd_initial_position.size())");


    double sum_dr_sq = 0.0;
    caviar2::Vector3d<int> domain_cross {0, 0, 0};
    int msd_particle_count = 0;
    for (unsigned int j = 0; j < pos_size; ++j)
    {
      if (msd_type > -1 && type[j] != static_cast<unsigned int>(msd_type))
        continue;

      auto dr = msd_initial_position[j] - pos[j];

      dr.x += atom_data->atom_struct_owned.msd_domain_cross[j].x * domain->size_global.x;
      dr.y += atom_data->atom_struct_owned.msd_domain_cross[j].y * domain->size_global.y;
      dr.z += atom_data->atom_struct_owned.msd_domain_cross[j].z * domain->size_global.z;

      domain_cross.x += atom_data->atom_struct_owned.msd_domain_cross[j].x;
      domain_cross.y += atom_data->atom_struct_owned.msd_domain_cross[j].x;
      domain_cross.z += atom_data->atom_struct_owned.msd_domain_cross[j].x;
      // domain-> periodic_distance(dr);

      auto dr_sq = dr * dr;
      sum_dr_sq += dr_sq;

      msd_particle_count++;
    }
    double msd = 0;
    if (msd_particle_count > 0)
      msd = sum_dr_sq / msd_particle_count;

    int domain_cross_tot = domain_cross.x + domain_cross.y+ domain_cross.z;
    if (my_mpi_rank == 0)
    {
      ofs_msd << i << " " << t << " " << msd << " " << domain_cross.x << " " << domain_cross.y << " " << domain_cross.z << " " << domain_cross_tot << std::endl;
    }
  }

} // writer

}
