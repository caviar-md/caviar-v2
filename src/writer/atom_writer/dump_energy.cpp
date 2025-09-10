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

#include "caviar2/domain.hpp"

namespace caviar2 {

namespace writer
{

  void Atom_writer::dump_energy_mpi(int64_t, double)
  {
    // double k_e = atom_data->kinetic_energy();

    // ofs_energy << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_writer::dump_energy_mpi_shared_atoms(int64_t i, double t)
  {
    double k_e = atom_data->kinetic_energy();

    if (my_mpi_rank != 0)
      return;

    ofs_energy << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_writer::dump_energy_serial(int64_t i, double t)
  {

#if defined(CAVIAR_WITH_MPI)
    if (energy_mpi_rank0)
#endif
    {
      double k_e = atom_data->kinetic_energy();

      ofs_energy << i << " " << t << " " << k_e << std::endl;
    }
  }
  void Atom_writer::dump_energy_mpi_per_process(int64_t i , double t)
  {

    if (energy_mpi_per_process)
    {
      double k_e = atom_data->kinetic_energy_mpi_domain();

      ofs_energy_mpi << i << " " << t << " " << k_e << std::endl;
    }

  }

  void Atom_writer::dump_temperature_mpi(int64_t, double)
  {
    // double k_e = atom_data->temperature();

    // ofs_temperature << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_writer::dump_temperature_mpi_shared_atoms(int64_t i, double t)
  {

    double k_e = atom_data->temperature();

    if (my_mpi_rank != 0)
      return;

    ofs_temperature << i << " " << t << " " << k_e << std::endl;
  }

  void Atom_writer::dump_temperature_serial(int64_t i, double t)
  {

#if defined(CAVIAR_WITH_MPI)
    if (temperature_mpi_rank0)
#endif
    {
      double k_e = atom_data->temperature();

      ofs_temperature << i << " " << t << " " << k_e << std::endl;
    }
  }
  void Atom_writer::dump_temperature_mpi_per_process(int64_t i, double t)
  {

    if (temperature_mpi_per_process)
    {
      double k_e = atom_data->temperature_mpi_domain();

      ofs_temperature_mpi << i << " " << t << " " << k_e << std::endl;
    }

  }

  void Atom_writer::dump_pressure_mpi(int64_t, double)
  {
    // double p = atom_data->pressure();

    // ofs_pressure << i << " " << t << " " << p << std::endl;
  }

  void Atom_writer::dump_pressure_mpi_shared_atoms(int64_t i, double t)
  {

    double p = atom_data->pressure();

    if (my_mpi_rank != 0)
      return;

    ofs_pressure << i << " " << t << " " << p << std::endl;
  }

  void Atom_writer::dump_pressure_serial(int64_t i, double t)
  {

#if defined(CAVIAR_WITH_MPI)
    if (pressure_mpi_rank0)
#endif
    {
      double p = atom_data->pressure();

      ofs_pressure << i << " " << t << " " << p << std::endl;
    }
  }

  void Atom_writer::dump_pressure_mpi_per_process(int64_t i, double t)
  {

    if (pressure_mpi_per_process)

    {
      double p = atom_data->pressure_mpi_domain();

      ofs_pressure_mpi << i << " " << t << " " << p << std::endl;
    }

  }

  void Atom_writer::dump_volume_mpi(int64_t, double)
  {
    // double p = atom_data->volume();

    // ofs_volume << i << " " << t << " " << p << std::endl;
  }

  void Atom_writer::dump_volume_mpi_shared_atoms(int64_t i, double t)
  {


    double v = domain->volume_global();

    if (my_mpi_rank != 0)
      return;

    ofs_volume << i << " " << t << " " << v << std::endl;
  }

  void Atom_writer::dump_volume_serial(int64_t i, double t)
  {

#if defined(CAVIAR_WITH_MPI)
    if (volume_mpi_rank0)
#endif
    {
      double v = domain->volume_global();

      ofs_volume << i << " " << t << " " << v << std::endl;
    }
  }

  void Atom_writer::dump_volume_mpi_per_process(int64_t i, double t)
  {

    if (volume_mpi_per_process)
    {
      double v = domain->volume_local();

      ofs_volume_mpi << i << " " << t << " " << v << std::endl;
    }

  }

} // writer

}
