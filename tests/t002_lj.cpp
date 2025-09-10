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

#include "caviar2/caviar2.hpp"
#include "caviar2/tools/all.hpp"
#include "caviar2/writer/atom_writer.hpp"

#include <iostream>

using namespace caviar2;

int main()
{
  Caviar2 caviar;

  // ===== global variables

  double cutoff_max = 5;

  double dt_all = 0.001;

  //===== Atom definition =====

  // #===== Domain

  caviar.domain.set_lower_global(Vector3d<double>{-51, -50, -50});
  caviar.domain.set_upper_global(Vector3d<double>{51, 50, 50});
  caviar.domain.set_boundary_condition(Vector3d<int>{1, 1, 1});
  caviar.domain.generate();

  // #========== Atom_data

  caviar.atom_data.set_ghost_cutoff(cutoff_max);

  for (int i = -20; i < 20; i += 5)
  {
    for (int j = -20; j < 20; j += 5)
    {
      for (int k = -20; k < 20; k += 5)
      {
        tools::Atom a1(&caviar);
        a1.set_type(0);
        a1.set_position(Vector3d<double>{(double)i, (double)j, (double)k});

        caviar.atom_data.add_atom(a1);
      }
    }
  }
  caviar.atom_data.add_type_mass(0, 1);
  caviar.atom_data.add_type_charge(0, 0);
  caviar.atom_data.set_k_b(1.0);
  caviar.atom_data.set_temperature_process(true);
  caviar.atom_data.add_random_velocity(1, .5);
  // #===== Neighborlist

  caviar.neighborlist.set_cutoff(cutoff_max);
  caviar.neighborlist.set_cutoff_extra_coef(1.0523);
  caviar.neighborlist.set_dt(dt_all);

  // #===== force_field

  force_field::Lj lj(&caviar);
  lj.set_cutoff(10);
  lj.set_sigma(0, 0, 1);
  lj.set_epsilon(0, 0, 1);

  caviar.add_force_field(&lj);

  // #====== writer
  writer::Atom_writer w1(&caviar);
  w1.set_xyz_step(500);
  w1.set_xyz_mpi_rank0(true);
  w1.set_xyz_mpi_per_process(true);

  w1.set_temperature_step(500);
  w1.set_temperature_mpi_rank0(false);
  w1.set_temperature_mpi_per_process(false);

  w1.set_energy_step(500);
  w1.set_energy_mpi_rank0(true);
  w1.set_energy_mpi_per_process(true);

  caviar.add_writer(&w1);

  // #=====  md simulator

  caviar.md_simulator.set_integrator_type(Integrator_t::Velocity_verlet);
  caviar.md_simulator.set_initial_step(0);
  caviar.md_simulator.set_final_step(200000);
  caviar.md_simulator.set_dt(dt_all);
  caviar.md_simulator.run();

  return 0;
}
