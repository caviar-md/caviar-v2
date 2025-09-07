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

#include "caviar2/forcefield.hpp"
#include <vector>

namespace caviar2
{

  namespace forcefield
  {

    /**
     * @brief Calculates Lennard-Jones (LJ) potential for particles.
     *
     * This class implements the Lennard-Jones force field calculations, including
     * options for lambda scaling, intra-molecular interaction filtering, and
     * support for WCA potential.
     */
    class Lj : public ForceField
    {
    public:
      /**
       * @brief Constructor.
       * @param c Pointer to the main CAVIAR object.
       */
      explicit Lj(class Caviar2 *caviar)
          : ForceField(caviar) // call base class constructor
      {
        // Your Lj-specific initialization here (if any)
      }

      /**
       * @brief Destructor.
       */
      ~Lj() {};

      /**
       * @brief Verifies that the settings are consistent and valid.
       *        Should be called after reading parameters.
       */
      void verify_settings();

      /**
       * @brief Calculates accelerations on particles due to LJ interactions.
       *        Overrides the pure virtual function in Force_field base class.
       */
      void calculate_acceleration();

    public:
      /// Lambda scaling factors for the attractive part of the potential [lambda_e(i,j)].
      std::vector<std::vector<double>> lambda_e;

      /// Lambda scaling factors for the repulsive part of the potential [lambda_s(i,j)].
      std::vector<std::vector<double>> lambda_s;

      /// Flag indicating whether lambda_e has been set by input.
      bool lambda_e_is_set = false;

      /// Flag indicating whether lambda_s has been set by input.
      bool lambda_s_is_set = false;

      /// If true, ignore intra-molecular interactions.
      bool ignore_intra_molecule = false;

      /// If true, input parameters are provided as arrays rather than per-atom.
      bool input_by_array;

      /// Epsilon parameters for each pair of types (depth of potential well).
      std::vector<std::vector<double>> epsilon;

      /// Sigma parameters for each pair of types (distance at which potential is zero).
      std::vector<std::vector<double>> sigma;

      /// If true, off-diagonal vectors (cross parameters) will be computed.
      bool make_off_diagonal_vectors;

      /// If true, input parameters are provided per atom type.
      bool input_by_atom;

      /// Epsilon parameters for individual atom types.
      std::vector<double> epsilon_atom;

      /// Sigma parameters for individual atom types.
      std::vector<double> sigma_atom;

      /// Debugging flag: enable jump fix algorithm.
      bool jump_fix;

      /// Debugging flag: monitor jump events.
      bool monitor_jump;

      /// Tolerance parameter for jump detection/fix.
      double jump_tol;

      /// If true, the Weeks-Chandler-Anderson (WCA) potential is activated.
      bool wca;

      /// If true, cutoff lists are activated for variable cutoff distances.
      bool cutoff_list_activated;

      /**
       * @brief List of cutoff distances between particle types.
       *        Used in specialized potentials like WCA.
       */
      std::vector<std::vector<double>> cutoff_list;
    };

  } // namespace force_field

} // namespace caviar2
