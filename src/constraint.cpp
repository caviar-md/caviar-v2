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

#include "caviar2/constraint.hpp"
#include "caviar2/caviar2.hpp"

namespace caviar2
{

  Constraint::Constraint(Caviar2 *fptr) : caviar_{fptr}
  {
    set_caviar(fptr);
  }

  Constraint::~Constraint()
  {
  }

  void Constraint::set_caviar(caviar2::Caviar2 *c)
  {
    caviar_ = c;
    atom_data = &(c->atom_data);
  }
  void Constraint::verify_settings()
  {
  }

  void Constraint::apply(int64_t) {}
  void Constraint::apply_shake(int64_t) {}
  void Constraint::fix_position(int64_t) {}
  void Constraint::fix_velocity(int64_t, bool &) {}
  void Constraint::apply_barostat(int64_t, bool &) {}
  void Constraint::apply_thermostat(int64_t, bool &) {}
  void Constraint::fix_acceleration(int64_t) {}

}
