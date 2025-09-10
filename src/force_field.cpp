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

#include "caviar2/force_field.hpp"
#include "caviar2/caviar2.hpp"

namespace caviar2
{

  Force_field::Force_field(Caviar2 *fptr) : caviar_{fptr}

  {
    set_caviar(caviar_);
  }
  
  void Force_field::set_caviar(Caviar2 *c) 
  { 
    caviar_ = c; 
    atom_data = &caviar_->atom_data;
    domain = &caviar_->domain;
    neighborlist = &caviar_->neighborlist;
  }

  void Force_field::verify_settings()
  {
  }

  double Force_field::energy()
  {
    caviar_->log.error_all(FC_FILE_LINE_FUNC, "The energy calculation of this force_field is not implemented");
    return 0.0;
  }

  double Force_field::potential(const Vector3d<double> &)
  {
    caviar_->log.error_all(FC_FILE_LINE_FUNC, "The potential calculation of this force_field is not implemented");
    return 0.0;
  }

  double Force_field::potential(const int)
  {
    caviar_->log.error_all(FC_FILE_LINE_FUNC, "The potential calculation of this force_field is not implemented");
    return 0.0;
  }

  Vector3d<double> Force_field::field(const Vector3d<double> &)
  {
    caviar_->log.error_all(FC_FILE_LINE_FUNC, "The field calculation of this force_field is not implemented");
    return Vector3d<double>{0, 0, 0};
  }

  Vector3d<double> Force_field::field(const int)
  {
    caviar_->log.error_all(FC_FILE_LINE_FUNC, "The field calculation of this force_field is not implemented");
    return Vector3d<double>{0, 0, 0};
  }

  void Force_field::scale_position(double, caviar2::Vector3d<int>)
  {
    caviar_->log.error_all(FC_FILE_LINE_FUNC, "The scale_position of this force_field is not implemented");
  }

}
