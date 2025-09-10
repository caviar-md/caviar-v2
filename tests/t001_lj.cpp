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
#include "caviar2/writer/atom_data.hpp"

#include <iostream>

using namespace caviar2;

int main()
{
  Caviar2 caviar;

  // ===== global variables

  double cutoff_max = 5;

  double dt_all = 0.001;

  // #===== Atom definition =====

  tools::Atom a1(&caviar);
  a1.type = 0;
  a1.position = Vector3d<double>{-1, 0, 0};

  tools::Atom a2(&caviar);
  ;
  a2.type = 0;
  a2.position = Vector3d<double>{1, 0, 0};

  // #===== Domain

  caviar.domain.lower_global = Vector3d<double>{-51, -50, -50};
  caviar.domain.upper_global = Vector3d<double>{51, 50, 50};
  caviar.domain.boundary_condition = Vector3d<int>{1, 1, 1};
  caviar.domain.generate();

  // #========== Atom_data

  caviar.atom_data.ghost_cutoff = cutoff_max;
  caviar.atom_data.add_atom(a1);
  caviar.atom_data.add_atom(a2);
  caviar.atom_data.add_type_mass(0, 1);
  caviar.atom_data.add_type_charge(0, 0);
  caviar.atom_data.set_k_b(1.0);
  caviar.atom_data.set_temperature_process(true);

  // #===== Neighborlist

  caviar.neighborlist.cutoff = cutoff_max;
  caviar.neighborlist.cutoff_extra_coef = 1.0523;
  caviar.neighborlist.dt = dt_all;

  // #===== force_field

  force_field::Lj lj(&caviar);
  lj.cutoff = 10;
  lj.set_sigma(0, 0, 1);
  lj.set_epsilon(0, 0, 1);

  caviar.add_force_field(&lj);

  // #====== writer
  writer::Atom_data w1(&caviar);
  w1.output_xyz = true;
  w1.xyz_step = 200;
  w1.xyz_mpi_rank0 = 1;
  w1.xyz_mpi_per_process = 1;
  w1.output_temperature = true;

  w1.temperature_step = 200;
  w1.temperature_mpi_rank0 = 0;
  w1.temperature_mpi_per_process = 0;
  w1.output_energy = true;

  w1.energy_step = 200;
  w1.energy_mpi_rank0 = 1;
  w1.energy_mpi_per_process = 1;

  caviar.add_writer(&w1);

  // #=====  md simulator

  caviar.md_simulator.set_integrator_type(Integrator_t::Velocity_verlet);
  caviar.md_simulator.set_initial_step(0);
  caviar.md_simulator.set_final_step(20000);
  caviar.md_simulator.set_dt(dt_all);
  caviar.md_simulator.run();

  return 0;
}
