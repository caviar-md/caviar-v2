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

namespace caviar2
{
  /**
   * More sharing results in less MPI overhead and probably faster MPI simulation. However, simulation will take more RAM.
   * For large number of processors and particles, it must be tested.
   */
  enum class MpiOptimization
  {
    None,
    SingleMdDomain,
    ShareAtoms,
    ShareMolecules,
    ShareAtomsAndMolecules,
  };
}