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

#include "caviar2/time_utility.hpp"
#include "caviar2/tools/time_function_3d.hpp"
#include "caviar2/domain.hpp"
#include "caviar2/communicator.hpp"
#include "caviar2/caviar2.hpp"

#include <ctime>
#include <sys/stat.h> // used for mkdir()

namespace caviar2
{

  namespace writer
  {

    Atom_writer::Atom_writer(Caviar2 *fptr) : Writer{fptr}
    {

      wallTimeXyzDump1 = get_wall_time();

      set_caviar(fptr);
    }

    Atom_writer::~Atom_writer()
    {
      if (ofs_xyz.is_open())
        ofs_xyz.close();

      if (ofs_xyz_mpi.is_open())
        ofs_xyz_mpi.close();

      if (ofs_xyz_ghost.is_open())
        ofs_xyz_ghost.close();

      if (ofs_xyz_ghost_mpi.is_open())
        ofs_xyz_ghost_mpi.close();

      if (ofs_energy.is_open())
        ofs_energy.close();

      if (ofs_energy_mpi.is_open())
        ofs_energy_mpi.close();

      if (ofs_temperature.is_open())
        ofs_temperature.close();

      if (ofs_temperature_mpi.is_open())
        ofs_temperature_mpi.close();

      if (ofs_pressure.is_open())
        ofs_pressure.close();

      if (ofs_pressure_mpi.is_open())
        ofs_pressure_mpi.close();

      if (ofs_povray.is_open())
        ofs_povray.close();

      if (ofs_povray_mpi.is_open())
        ofs_povray_mpi.close();

      if (ofs_msd.is_open())
        ofs_msd.close();

      if (ofs_msd_mpi.is_open())
        ofs_msd_mpi.close();

      if (ofs_volume.is_open())
        ofs_volume.close();

      if (ofs_volume_mpi.is_open())
        ofs_volume_mpi.close();
    }

    void Atom_writer::set_caviar(caviar2::Caviar2 *c)
    {
      caviar_ = c;
      atom_data = &(c->atom_data);
    }


   
    void Atom_writer::initialize()
    {
      initialized = true;

      if (output_msd)
      {
        if (!atom_data->get_msd_process())
          caviar_->log.error_all(FC_FILE_LINE_FUNC, "In order to have 'output_msd' in writer::Atom_data, 'msd_process' must be activated in atom_data::Atom_data");
      }

      if (output_volume)
      {
      }
      // --- just to make povray outpuy folder ---
      if (my_mpi_rank == 0)
      {
        if (output_povray)
        {
          std::string str_folder_pov;
          str_folder_pov.append("o_pov");
          const char *char_folder_pov = str_folder_pov.c_str();
          mkdir(char_folder_pov, 0777); // make povray  output folder //
        }
      }
      /*
      #ifdef CAVIAR_WITH_MPI
        sprintf ( buffer, "_me%u", comm->me );
        str_filename.append ( buffer);
      #endif
      */

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
      if (my_mpi_rank == 0)
      {
        if (output_xyz)
          if (!ofs_xyz.is_open())
            ofs_xyz.open((file_name_xyz + ".xyz").c_str());

        if (output_xyz_ghost)
          if (!ofs_xyz_ghost.is_open())
            ofs_xyz_ghost.open((file_name_xyz_ghost + ".xyz").c_str());

        if (output_povray)
          if (!ofs_povray.is_open())
            ofs_povray.open((file_name_povray + ".pov").c_str());

        if (output_energy)
          if (!ofs_energy.is_open())
            ofs_energy.open((file_name_energy + ".txt").c_str());

        if (output_temperature)
          if (!ofs_temperature.is_open())
            ofs_temperature.open((file_name_temperature + ".txt").c_str());

        if (output_pressure)
          if (!ofs_pressure.is_open())
            ofs_pressure.open((file_name_pressure + ".txt").c_str());

        if (output_msd)
          if (!ofs_msd.is_open())
            ofs_msd.open((file_name_msd + ".txt").c_str());

        if (output_volume)
          if (!ofs_volume.is_open())
            ofs_volume.open((file_name_volume + ".txt").c_str());
      }
#elif defined(CAVIAR_WITH_MPI)

      if (output_xyz && xyz_mpi_per_process)
        if (!ofs_xyz_mpi.is_open())
          ofs_xyz_mpi.open((file_name_xyz + "_mpi" + std::to_string(my_mpi_rank) + ".xyz").c_str());

      if (output_xyz_ghost && xyz_ghost_mpi_per_process)
        if (!ofs_xyz_ghost_mpi.is_open())
          ofs_xyz_ghost_mpi.open((file_name_xyz_ghost + "_mpi" + std::to_string(my_mpi_rank) + ".xyz").c_str());

      if (output_povray && povray_mpi_per_process)
        if (!ofs_povray_mpi.is_open())
          ofs_povray_mpi.open((file_name_povray + "_mpi" + std::to_string(my_mpi_rank) + ".pov").c_str());

      if (output_energy && energy_mpi_per_process)
        if (!ofs_energy_mpi.is_open())
          ofs_energy_mpi.open((file_name_energy + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

      if (output_temperature && temperature_mpi_per_process)
        if (!ofs_temperature_mpi.is_open())
          ofs_temperature_mpi.open((file_name_temperature + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

      if (output_pressure && pressure_mpi_per_process)
        if (!ofs_pressure_mpi.is_open())
          ofs_pressure_mpi.open((file_name_pressure + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

      if (output_msd && msd_mpi_per_process)
        if (!ofs_msd_mpi.is_open())
          ofs_msd_mpi.open((file_name_msd + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

      if (output_volume && volume_mpi_per_process)
        if (!ofs_volume_mpi.is_open())
          ofs_volume_mpi.open((file_name_volume + "_mpi" + std::to_string(my_mpi_rank) + ".txt").c_str());

      if (my_mpi_rank == 0)
      {

        if (output_xyz)
          if (!ofs_xyz.is_open())
            ofs_xyz.open((file_name_xyz + ".xyz").c_str());

        if (output_xyz_ghost)
          if (!ofs_xyz_ghost.is_open())
            ofs_xyz_ghost.open((file_name_xyz_ghost + ".xyz").c_str());

        if (output_povray)
          if (!ofs_povray.is_open())
            ofs_povray.open((file_name_povray + ".pov").c_str());

        if (output_energy)
          if (!ofs_energy.is_open())
            ofs_energy.open((file_name_energy + ".txt").c_str());

        if (output_temperature)
          if (!ofs_temperature.is_open())
            ofs_temperature.open((file_name_temperature + ".txt").c_str());

        if (output_pressure)
          if (!ofs_pressure.is_open())
            ofs_pressure.open((file_name_pressure + ".txt").c_str());

        if (output_msd)
          if (!ofs_msd.is_open())
            ofs_msd.open((file_name_msd + ".txt").c_str());

        if (output_volume)
          if (!ofs_volume.is_open())
            ofs_volume.open((file_name_volume + ".txt").c_str());
      }
#else

      if (output_xyz)
        if (!ofs_xyz.is_open())
          ofs_xyz.open((file_name_xyz + ".xyz").c_str());

      if (output_xyz_ghost)
        if (!ofs_xyz_ghost.is_open())
          ofs_xyz_ghost.open((file_name_xyz_ghost + ".xyz").c_str());

      if (output_povray)
        if (!ofs_povray.is_open())
          ofs_povray.open((file_name_povray + ".pov").c_str());

      if (output_energy)
        if (!ofs_energy.is_open())
          ofs_energy.open((file_name_energy + ".txt").c_str());

      if (output_temperature)
        if (!ofs_temperature.is_open())
          ofs_temperature.open((file_name_temperature + ".txt").c_str());

      if (output_pressure)
        if (!ofs_pressure.is_open())
          ofs_pressure.open((file_name_pressure + ".txt").c_str());

      if (output_msd)
        if (!ofs_msd.is_open())
          ofs_msd.open((file_name_msd + ".txt").c_str());

      if (output_volume)
        if (!ofs_volume.is_open())
          ofs_volume.open((file_name_volume + ".txt").c_str());
#endif
    }

    void Atom_writer::open_files() {}
    void Atom_writer::close_files() {}
    void Atom_writer::generate() {}

    void Atom_writer::report_xyz_dump(int64_t i, double)
    {

      // if (my_mpi_rank != 0)
      //   return;

      double wallTimeXyzDump2 = get_wall_time();

      double dtstart = wallTimeXyzDump2 - wallTimeXyzDump1;
      std::string s = "dump_xyz [" + std::to_string(i) +
                      +"] (" + std::to_string(dtstart) + " S)";

      caviar_->log.info(s, 2);

      wallTimeXyzDump1 = wallTimeXyzDump2;
    }

    void Atom_writer::write(int64_t i, double t)
    {

      if (!initialized)
        initialize();

      report_xyz_this_time_step = false;

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

      if (my_mpi_rank == 0)
      {
        write_serial(i, t);
      }
#elif defined(CAVIAR_WITH_MPI)

      write_mpi_shared_atoms(i, t);

      write_mpi_per_process(i, t);

#else

      write_serial(i, t);

#endif
      if (report_xyz_this_time_step)
        report_xyz_dump(i, t);
    }

    void Atom_writer::write_mpi(int64_t i, double t)
    {
      if (output_xyz)
        if (i % xyz_step == 0)
        {
          dump_xyz_mpi(i, t);
          report_xyz_this_time_step = true;
        }

      if (output_xyz_ghost)
        if (i % xyz_ghost_step == 0)
          dump_xyz_ghost_mpi(i, t);

      if (output_energy)
        if (i % energy_step == 0)
          dump_energy_mpi(i, t);

      if (output_temperature)
        if (i % temperature_step == 0)
          dump_temperature_mpi(i, t);

      if (output_pressure)
        if (i % pressure_step == 0)
          dump_pressure_mpi(i, t);

      if (output_povray)
        if (i % povray_step == 0)
          dump_povray_mpi(i, t);

      if (output_msd)
        if (i % msd_step == 0)
          dump_msd_mpi(i, t);

      if (output_volume)
        if (i % volume_step == 0)
          dump_volume_mpi(i, t);
    }

    // All MPI processes work on this but only rank0 outputs to a file
    void Atom_writer::write_mpi_shared_atoms(int64_t i, double t)
    {
      if (output_xyz && xyz_mpi_rank0)
        if (i % xyz_step == 0)
        {
          dump_xyz_mpi_shared_atoms(i, t);
          report_xyz_this_time_step = true;
        }

      if (output_xyz_ghost && xyz_ghost_mpi_rank0)
        if (i % xyz_ghost_step == 0)
          dump_xyz_ghost_mpi_shared_atoms(i, t);

      if (output_energy && energy_mpi_rank0)
        if (i % energy_step == 0)
          dump_energy_mpi_shared_atoms(i, t);

      if (output_temperature && temperature_mpi_rank0)
        if (i % temperature_step == 0)
          dump_temperature_mpi_shared_atoms(i, t);

      if (output_pressure && pressure_mpi_rank0)
        if (i % pressure_step == 0)
          dump_pressure_mpi_shared_atoms(i, t);

      if (output_povray && povray_mpi_rank0)
        if (i % povray_step == 0)
          dump_povray_mpi_shared_atoms(i, t);

      if (output_msd && msd_mpi_rank0)
        if (i % msd_step == 0)
          dump_msd_mpi_shared_atoms(i, t);

      if (output_volume && volume_mpi_rank0)
        if (i % volume_step == 0)
          dump_volume_mpi_shared_atoms(i, t);
    }

    // All MPI processes work on this but only rank0 outputs to a file
    void Atom_writer::write_mpi_per_process(int64_t i, double t)
    {
      if (output_xyz && xyz_mpi_per_process)
        if (i % xyz_step == 0)
        {
          dump_xyz_mpi_per_process(i, t);
          report_xyz_this_time_step = true;
        }

      if (output_xyz_ghost && xyz_ghost_mpi_per_process)
        if (i % xyz_ghost_step == 0)
          dump_xyz_ghost_mpi_per_process(i, t);

      if (output_energy && energy_mpi_rank0)
        if (i % energy_step == 0)
          dump_energy_mpi_per_process(i, t);

      if (output_temperature && temperature_mpi_per_process)
        if (i % temperature_step == 0)
          dump_temperature_mpi_per_process(i, t);

      if (output_pressure && pressure_mpi_per_process)
        if (i % pressure_step == 0)
          dump_pressure_mpi_per_process(i, t);

      if (output_povray && povray_mpi_per_process)
        if (i % povray_step == 0)
          dump_povray_mpi_per_process(i, t);

      if (output_msd && msd_mpi_per_process)
        if (i % msd_step == 0)
          dump_msd_mpi_per_process(i, t);

      if (output_volume && volume_mpi_per_process)
        if (i % volume_step == 0)
          dump_volume_mpi_per_process(i, t);
    }

    void Atom_writer::write_serial(int64_t i, double t)
    {
      if (output_xyz)
        if (i % xyz_step == 0)
        {
          dump_xyz_serial(i, t);
          report_xyz_this_time_step = true;
        }

      if (output_xyz_ghost)
        if (i % xyz_ghost_step == 0)
          dump_xyz_ghost_serial(i, t);

      if (output_energy)
        if (i % energy_step == 0)
          dump_energy_serial(i, t);

      if (output_temperature)
        if (i % temperature_step == 0)
          dump_temperature_serial(i, t);

      if (output_pressure)
        if (i % pressure_step == 0)
          dump_pressure_serial(i, t);

      if (output_povray)
        if (i % povray_step == 0)
          dump_povray_serial(i, t);

      if (output_msd)
        if (i % msd_step == 0)
          dump_msd_serial(i, t);

      if (output_volume)
        if (i % volume_step == 0)
          dump_volume_serial(i, t);
    }

    void Atom_writer::start_new_files() {}              // add_time_to_previous
    void Atom_writer::start_new_files(std::string &) {} // add_time_to_previous

  } // Atom_data

}
