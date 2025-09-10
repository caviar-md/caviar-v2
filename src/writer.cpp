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

#include "caviar2/writer.hpp"
#include "caviar2/caviar2.hpp"

namespace caviar2
{

    Writer::Writer(Caviar2 *fptr) : caviar_{fptr}, initialized{false}, my_mpi_rank{fptr->comm.me},
                                    mpi_world_size{fptr->comm.nprocs}
    {
    }

    Writer::~Writer()
    {
    }

    void Writer::set_caviar(caviar2::Caviar2 *c)
    {
        caviar_ = c;
    }

    void Writer::initialize() {}
    void Writer::write(int64_t, double) {}         // time_step and time
    void Writer::start_new_files() {}              // add_time_to_previous
    void Writer::start_new_files(std::string &) {} // add_time_to_previous
    void Writer::open_files() {}
    void Writer::close_files() {}
    void Writer::generate() {}

}
