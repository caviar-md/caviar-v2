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

#include "caviar2/writer/force_writer.hpp"

// #include "caviar2/atom_data.hpp"
// #include <ctime>

namespace caviar2
{

  namespace writer
  {

    Force_writer::Force_writer(Caviar2 *fptr) : Writer{fptr}
    {
    }

    Force_writer::~Force_writer()
    {
    }
    /*
      bool Force_writer::read(caviar2::interpreter::Parser *parser)
      {
        FC_OBJECT_READ_INFO
        bool in_file = true;
        while (true)
        {
          GET_A_TOKEN_FOR_CREATION
          auto t = token.string_value;
          FC_OBJECT_READ_INFO_STR
          if (string_cmp(t, "output_all_acc"))
          {
            write();
          }  else  if (string_cmp(t,"set_atom_data") || string_cmp(t,"atom_data")) {
             FIND_OBJECT_BY_NAME(atom_data,it)
             atom_data = object_container->atom_data[it->second.index];
           }
        }
        return in_file;
      }
    */
    void Force_writer::initialize() {}
    void Force_writer::write()
    {
      /*std::ofstream ofs ("o_acc");
      const auto &pos = atom_data -> atom_struct_owned.position;
      const auto &acc = atom_data -> atom_struct_owned.acceleration;
      for (unsigned int i=0;i<pos.size();++i) {
        ofs << i << " " << acc[i].x << "\t" << acc[i].y << "\t" << acc[i].z << "\n" ;
      }*/
    }
    void Force_writer::write(int64_t) {}                 // current time_step
    void Force_writer::write(double) {}                  // current time
    void Force_writer::write(int64_t, double) {}         // time_step and time
    void Force_writer::start_new_files() {}              // add_time_to_previous
    void Force_writer::start_new_files(std::string &) {} // add_time_to_previous

  } // Force_field

}
