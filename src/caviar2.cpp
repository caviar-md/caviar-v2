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

namespace caviar2
{

    Caviar2::Caviar2()
    {
        atom_data.set_caviar(this);
        domain.set_caviar(this);
        md_simulator.set_caviar(this);
        neighborlist.set_caviar(this);
        log.set_caviar(this);
        comm.set_caviar(this);
    }

    Caviar2::~Caviar2()
    {
        // for (auto f : forces)
        //     delete f;
        // for (auto c : constraints)
        //     delete c;
        // for (auto w : writers)
        //     delete w;
    }

    void Caviar2::add_force_field(Force_field *x)
    {
        log.info("Force_field added to the list");
        forces.push_back(x);
    }

    void Caviar2::add_constraint(Constraint *x)
    {
        constraints.push_back(x);
    }

    void Caviar2::add_writer(Writer *x)
    {
        writers.push_back(x);
    }

} // namespace caviar2
