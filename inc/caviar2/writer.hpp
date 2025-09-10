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

namespace caviar2
{

  /**
   * This class is the base class for all the writers.
   *
   *
   */
  class Writer
  {
  public:
    Writer();
    Writer(class Caviar2* caviar);
    virtual ~Writer();
    virtual void initialize();
    virtual void write(int64_t, double);         // time_step and time
    virtual void start_new_files();              // add_time_to_previous
    virtual void start_new_files(std::string &); // add_time_to_previous
    virtual void open_files();
    virtual void close_files();
    virtual void generate();
    bool initialized;
    int64_t last_timestep;
    double last_time;
    double dt;
    int my_mpi_rank, mpi_world_size;
    // FC_BASE_OBJECT_COMMON_TOOLS
    virtual void set_caviar(Caviar2 *c);

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

}