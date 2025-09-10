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

#include "caviar2/writer.hpp"

namespace caviar2 {

namespace writer
{

    /**
     * This class has a writer for force-field
     *
     *
     */
    class Force_writer : public Writer
    {
    public:
        Force_writer(class Caviar2* caviar);
        Force_writer();
        ~Force_writer();
        
        void initialize();
        void write();
        void write(int64_t);                 // current time_step
        void write(double);                  // current time
        void write(int64_t, double);         // time_step and time
        void start_new_files();              // add_time_to_previous
        void start_new_files(std::string &); // add_time_to_previous

    public:
       virtual void set_caviar(Caviar2 *c);

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership

    };

} // writer

}

