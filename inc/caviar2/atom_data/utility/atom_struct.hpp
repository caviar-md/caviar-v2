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

#include <vector>
#include "caviar2/vector3d.hpp"

namespace caviar2
{

  namespace atom_data
  {
    /**
     * It contains all the physical data of atoms.
     */
    struct Atom_struct
    {
      /**
       * MPI rank of the atom. In this case, all processes have a copy of the atom
       */
      std::vector<size_t> mpi_rank;

      /**
       * 'id' is a global and tools number assigned to an atom.
       */
      std::vector<size_t> id;

      // /**
      //  * A tag to be set on the atoms.
      //  */
      // std::vector<size_t> tag;

      /**
       * Atom type decides the charge, mass and any other property shared between
       * a defined type (for example, Elements).
       */
      std::vector<size_t> type;

      /**
       * Different by atom_id. Can be changed in simulation (it is needed to be sent-recv. by MPI)
       */
      std::vector<double> charge_atom;

      /**
       * Different by atom_id. Can be changed in simulation (it is needed to be sent-recv. by MPI)
       */
      std::vector<double> mass_atom;

      /**
       * Atom kinematic properties in the current time-step.
       */
      std::vector<Vector3d<double>> position, velocity, acceleration;

      /**
       * This vectors are used in some integrator schemes and constraint methods.
       * They can be defined in their related objects, but they may be needed in
       * more than one objects at once (for example, constraint::M_shake and
       * integrator::Leap_frog). This makes it the reason to define it here.
       * This function may be needed to have MPI_send-recv. process in these case.
       * look up to it.
       */
      std::vector<Vector3d<double>> position_old, velocity_old, acceleration_old;

      /**
       * Coordinates before applying any constraints. It can be used to calculate constraint forces which is needed
       * in pressure calculations.
       */
      // std::vector<Vector3d<double>> position_no_constraint, velocity_no_constraint, acceleration_no_constraint;

      /**
       * this vector is meaningful when there's one domain. We can calculate MSD
       * using this. It collects number of periodic domain cross for each particle.
       */
      std::vector<Vector3d<int32_t>> msd_domain_cross;

      /**
       * this vector contain a molecule index for all the atoms. if it's '-1' the
       * atom is not of any molecule. This matters in the MPI process. All of the
       * atoms of a molecule should be existed in one process.
       */
      std::vector<int32_t> molecule_index; //

      /**
       * Number of atomic bonds each atom have. It is used to limit
       * the bond creations.
       */
      std::vector<size_t> atomic_bond_count;
    };
  }
}