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

#include <string>

#include "caviar2/atom_data.hpp"
#include "caviar2/constraint.hpp"
#include "caviar2/domain.hpp"
#include "caviar2/neighborlist.hpp"
#include "caviar2/force_field.hpp"
#include "caviar2/md_simulator.hpp"
#include "caviar2/writer.hpp"
#include "caviar2/log.hpp"
#include "caviar2/communicator.hpp"

#include "caviar2/force/lj.hpp"

namespace caviar2
{

  /**
   * @brief A sample class for caviar2 library.
   *
   * This class demonstrates a basic interface.
   */
  class Caviar2
  {
  public:
    Caviar2();
    ~Caviar2();


    Atom_data atom_data;
    Domain domain;
    Md_simulator md_simulator;
    Neighborlist neighborlist;
    Log log;
    Communicator comm;

    std::vector<Force_field *> forces;
    std::vector<Constraint *> constraints;
    std::vector<Writer *> writers;

    void add_force_field(Force_field *x);

    void add_constraint(Constraint *x);

    void add_writer(Writer *x);

  private:

    // Private implementation details (PIMPL or just members)
  };

} // namespace caviar2
