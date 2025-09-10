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

#include <inttypes.h>

namespace caviar2
{

  class Integrator;
  class Atom_data;

  enum class Constraint_t
  {
    Atom_molarity,
    Atoms_molarity,
    Berendsen,
    Cm_motion,
    M_shake,
    Nose_hoover,
    Nve,
    Rattle,
    Shake,
    Unknown
  };

  /**
   * This class is the base class for all the constraints.
   * A constraint is called during a simulations. It can be implemented in
   * different parts of the virtual functions to be called at the desired moments
   * see the md_simulator objects to understant the time of call.
   */
  class Constraint
  {
  public:
    Constraint();
    Constraint(class Caviar2 *caviar);
    virtual ~Constraint();

    /**
     *   it will be applied at the end of each time step
     */
    virtual void apply(int64_t);

    /**
     * it should be used after calculating new position to fix it before using it for acceleration calculation.
     * for example, shake, m-shake, rattle use it
     */
    virtual void fix_position(int64_t);

    /**
     * it should be used after calculating velocity to fix it.
     */
    virtual void fix_velocity(int64_t, bool &recalculate_temperature);

    /**
     * it should be used after calculating acceleration to fix it.
     */
    virtual void fix_acceleration(int64_t);

    /**
     * shake, m-shake, rattle use it
     */
    virtual void apply_shake(int64_t);

    /**
     * NVE, NVT fix. It applies on velocity, but it must be called only on specific algorithm step
     */
    virtual void apply_thermostat(int64_t, bool &recalculate_temperature);

    /**
     * NPT fix. It applies on position, but it must be called only on specific algorithm step.
     * if (fix_position_needed) returns true, it means that the barostat scaling is applied
     * on position and one needs to call shake algorithm again
     */
    virtual void apply_barostat(int64_t, bool &fix_position_needed);

    class Atom_data *atom_data = nullptr;

    Constraint_t constraint_type;

    /**
     * MPI rank of the classs
     */
    int my_mpi_rank = -1;

    void verify_settings();

  public:
    virtual void set_caviar(Caviar2 *c) ;
  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

}