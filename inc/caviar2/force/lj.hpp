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

#include "caviar2/force_field.hpp"
#include <vector>

namespace caviar2
{

  namespace force_field
  {

    /**
     * @brief Calculates Lennard-Jones (LJ) potential for particles.
     *
     * This class implements the Lennard-Jones force field calculations, including
     * options for lambda scaling, intra-molecular interaction filtering, and
     * support for WCA potential.
     */
    class Lj : public Force_field
    {
    public:
      /**
       * @brief Constructor.
       * @param c Pointer to the main CAVIAR object.
       */
      Lj(class Caviar2 *caviar);

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




    // Single-value parameters
    void set_cutoff(double value);
    void set_jump_tol(double value);

    // Flags
    void enable_monitor_jump(bool flag = true);
    void ignore_intra_molecule_flag(bool flag = true);
    void enable_wca(bool flag = true);
    void enable_off_diagonal_vectors(bool flag = true);

    // Array / vector parameters
    void set_cutoff_list(const std::vector<std::vector<double>>& values);
    void set_epsilon(const int type1, const int type2, const double value);
    void set_sigma(const int type1, const int type2, const double value);
    void set_lambda_s(const std::vector<std::vector<double>>& values);
    void set_lambda_e(const std::vector<std::vector<double>>& values);
    void set_epsilon_atom(const int type, const double value);
    void set_sigma_atom(const int type, const double value);

    // References to other objects
    void set_neighborlist(class NeighborList* nl);
    void set_atom_data(class Atom_data* ad);
    void set_domain(class Domain* dom);
    
    private:
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
      bool input_by_array = false;

      /// Epsilon parameters for each pair of types (depth of potential well).
      std::vector<std::vector<double>> epsilon;

      /// Sigma parameters for each pair of types (distance at which potential is zero).
      std::vector<std::vector<double>> sigma;

      /// If true, off-diagonal vectors (cross parameters) will be computed.
      bool make_off_diagonal_vectors = false;

      /// If true, input parameters are provided per atom type.
      bool input_by_atom = false;

      /// Epsilon parameters for individual atom types.
      std::vector<double> epsilon_atom;

      /// Sigma parameters for individual atom types.
      std::vector<double> sigma_atom;

      /// Debugging flag: enable jump fix algorithm.
      bool jump_fix = false;

      /// Debugging flag: monitor jump events.
      bool monitor_jump = false;

      /// Tolerance parameter for jump detection/fix.
      double jump_tol = 1e-6;

      /// If true, the Weeks-Chandler-Anderson (WCA) potential is activated.
      bool wca = false;

      /// If true, cutoff lists are activated for variable cutoff distances.
      bool cutoff_list_activated = false;

      /**
       * @brief List of cutoff distances between particle types.
       *        Used in specialized potentials like WCA.
       */
      std::vector<std::vector<double>> cutoff_list;
    };

  } // namespace force_field

} // namespace caviar2
