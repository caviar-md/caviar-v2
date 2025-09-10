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

#include "caviar2/shape.hpp"

namespace caviar2
{

  Shape::Shape(Caviar2 *fptr)
  {
  }

  Shape::~Shape()
  {
  }

  void Shape::verify_settings()
  {
  }

  bool Shape::is_outside(const Vector3d<double> &v)
  {
    return !is_inside(v);
  }

  bool Shape::is_outside(const Vector3d<double> &v, const double r)
  {
    return !is_inside(v, r);
  }

  void Shape::scale_position(double, Vector3d<int>)
  {
    // error->all(FC_FILE_LINE_FUNC, "The scale_position of this force_field is not implemented");
  }
}
