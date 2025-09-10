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
#include "caviar2/vector3d.hpp"

namespace caviar2
{
  class Caviar2;

  inline void normalize(Vector3d<double> &v)
  {
    v /= std::sqrt(v * v);
  }

  /**
   * This class is the base class for all the shapes.
   *
   *
   */
  class Shape
  {
  public:
    /**
     * Constructor.
     */
    Shape(class Caviar2 *);
    Shape();

    /**
     * Destructor.
     */
    virtual ~Shape();

    virtual bool is_inside(const Vector3d<double> &) = 0;
    virtual bool is_outside(const Vector3d<double> &);
    virtual bool is_inside(const Vector3d<double> &, const double rad) = 0;
    virtual bool is_outside(const Vector3d<double> &, const double rad);
    virtual bool in_contact(const Vector3d<double> &, const double rad, Vector3d<double> &contact_vector) = 0;

    /**
     * Used in barostat scaling for geometrical forces.
     *
     */
    virtual void scale_position(double scale_ratio, caviar2::Vector3d<int> scale_axis);
    virtual void set_caviar(Caviar2 *c);
    void verify_settings();

  protected:
    Caviar2 *caviar_;
  };
}