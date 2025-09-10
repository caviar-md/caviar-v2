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

#include "caviar2/vector3d.hpp"
#include "caviar2/writer.hpp"
#include "caviar2/atom_data/utility/mpi_packet_info.hpp"
#include <fstream>
#include <vector>

namespace caviar2 {

class Atom_data;
class Domain;
namespace tools
{
  class Time_function_3d;
}
namespace writer
{

  /**
   * This class has a writer for Atom_data class
   *
   *
   */
  class Atom_writer : public Writer
  {
  public:
    Atom_writer(class Caviar2* caviar) ;
    Atom_writer();
    ~Atom_writer();
    
    void initialize();

    void write(int64_t, double);         // time_step and time
    void write_mpi(int64_t, double);         // time_step and time
    void write_mpi_per_process(int64_t, double);         // time_step and time
    void write_mpi_shared_atoms(int64_t, double);         // time_step and time
    void write_serial(int64_t, double);  // time_step and time

    void start_new_files();              // add_time_to_previous
    void start_new_files(std::string &); // add_time_to_previous
    void open_files();
    void close_files();
    void generate();
    void dump_kinetic_energy();

    tools::Time_function_3d *position_offset = nullptr;

    /////////////

    caviar2::Atom_data *atom_data = nullptr;

    // used in msd calculations
    Domain *domain = nullptr;

    void dump_energy_mpi(int64_t, double); // dump energy to file
    void dump_energy_mpi_per_process(int64_t, double); // dump energy to file
    void dump_energy_mpi_shared_atoms(int64_t, double); // dump energy to file
    void dump_energy_serial(int64_t, double); // dump energy to file

    void dump_temperature_mpi(int64_t, double); // dump temperature to file
    void dump_temperature_mpi_per_process(int64_t, double); // dump temperature to file
    void dump_temperature_mpi_shared_atoms(int64_t, double); // dump temperature to file
    void dump_temperature_serial(int64_t, double); // dump temperature to file

    void dump_pressure_mpi(int64_t, double); // dump pressure to file
    void dump_pressure_mpi_per_process(int64_t, double); // dump pressure to file
    void dump_pressure_mpi_shared_atoms(int64_t, double); // dump pressure to file
    void dump_pressure_serial(int64_t, double); // dump pressure to file

    void dump_xyz_mpi(int64_t, double);             // dump positions to file in xyz format. Number of atoms are different in mpi domains
    void dump_xyz_mpi_per_process(int64_t, double);             // dump positions to file in xyz format. Number of atoms are different in mpi domains
    void dump_xyz_mpi_shared_atoms(int64_t, double);  // dump positions to file in xyz format. Number of atoms are the same in mpi domains
    void dump_xyz_serial(int64_t, double);   // dump positions to file in xyz format.

    void dump_xyz_ghost_mpi(int64_t, double );        // dump positions to file in xyz format. Number of atoms are different in mpi domains
    void dump_xyz_ghost_mpi_per_process(int64_t, double );        // dump positions to file in xyz format. Number of atoms are different in mpi domains
    void dump_xyz_ghost_mpi_shared_atoms(int64_t, double );    // dump positions to file in xyz format. Number of atoms are different in mpi domains
    void dump_xyz_ghost_serial(int64_t, double);      // dump positions to file in xyz format

    void dump_povray_mpi(int64_t, double); // dump positions to snapshot files in povray format
    void dump_povray_mpi_per_process(int64_t, double); // dump positions to snapshot files in povray format
    void dump_povray_mpi_shared_atoms(int64_t, double); // dump positions to snapshot files in povray format
    void dump_povray_serial(int64_t, double); // dump positions to snapshot files in povray format

    void dump_msd_mpi(int64_t, double); //
    void dump_msd_mpi_per_process(int64_t, double); //
    void dump_msd_mpi_shared_atoms(int64_t, double); //
    void dump_msd_serial(int64_t, double); //


    void dump_volume_mpi(int64_t, double); // dump volume to file
    void dump_volume_mpi_per_process(int64_t, double); // dump volume to file
    void dump_volume_mpi_shared_atoms(int64_t, double); // dump volume to file
    void dump_volume_serial(int64_t, double); // dump volume to file

    void report_xyz_dump(int64_t, double);

    bool report_xyz_this_time_step = false;
    
    int64_t energy_step = 100;
    int64_t temperature_step = 100;
    int64_t pressure_step = 100;
    int64_t xyz_step = 100;
    int64_t xyz_ghost_step = 100;
    int64_t povray_step = 100;
    int64_t msd_step = 100; // number of steps to output data
    int64_t volume_step = 100;

    int msd_type = -1; // type of atom that should be used in msd calculations. If '-1' all of the atoms will be calculated
    int msd_initial_step = 0; // At this step, msd_initial_position will be generated

    std::ofstream ofs_energy;
    std::ofstream ofs_energy_mpi;
    std::ofstream ofs_temperature;
    std::ofstream ofs_temperature_mpi;
    std::ofstream ofs_pressure;
    std::ofstream ofs_pressure_mpi;
    std::ofstream ofs_xyz;
    std::ofstream ofs_xyz_mpi;
    std::ofstream ofs_xyz_ghost;
    std::ofstream ofs_xyz_ghost_mpi;
    std::ofstream ofs_povray;
    std::ofstream ofs_povray_mpi;
    std::ofstream ofs_msd; // mean square distance
    std::ofstream ofs_msd_mpi; // mean square distance
    std::ofstream ofs_volume; // mean square distance
    std::ofstream ofs_volume_mpi; // mean square distance

    // if true, outputs would be created    
    bool output_energy = false;
    bool output_temperature = false;
    bool output_pressure = false;
    bool output_xyz = false;
    bool output_xyz_ghost = false;
    bool output_povray = false;
    bool output_msd = false;
    bool output_volume = false;
    
    std::string file_name_xyz = "o_xyz";
    std::string file_name_xyz_ghost = "o_xyz_g";
    std::string file_name_energy = "o_energy";
    std::string file_name_temperature = "o_temperature";
    std::string file_name_pressure = "o_pressure";
    std::string file_name_povray = "o_pov";
    std::string file_name_msd = "o_msd";
    std::string file_name_volume = "o_volume";

    // dump velocity and acceleration alongside position in xyz file
    bool xyz_output_velocity = false; 
    bool xyz_output_acceleration = false; 
    bool xyz_output_id = false; 

    bool xyz_mpi_per_process = false;
    bool xyz_ghost_mpi_per_process = false;
    bool energy_mpi_per_process = false;
    bool temperature_mpi_per_process = false;
    bool pressure_mpi_per_process = false;
    bool povray_mpi_per_process = false;
    bool msd_mpi_per_process = false;
    bool volume_mpi_per_process = false;

    bool xyz_mpi_rank0 = false;
    bool xyz_ghost_mpi_rank0 = false;
    bool energy_mpi_rank0 = false;
    bool temperature_mpi_rank0 = false;
    bool pressure_mpi_rank0 = false;
    bool povray_mpi_rank0 = false;
    bool msd_mpi_rank0 = false;
    bool volume_mpi_rank0 = false;

    std::vector<caviar2::Vector3d<double>> msd_initial_position;

    // records previous wallTime of XYZ dump.
    double wallTimeXyzDump1;

    /////////

  /**
   * It stores the info about location of the owned Atom_struct data inside the MPI packets.
   */
  atom_data::MPI_packet_info mpi_packet_info_owned;

  /**
   * It stores the info about location of the ghost Atom_struct data inside the MPI packets.
   */
  atom_data::MPI_packet_info mpi_packet_info_ghost;
  

  public:
     void set_caviar(Caviar2 *c);

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // writer

}

