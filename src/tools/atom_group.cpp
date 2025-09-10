
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// This file is part of the CAVIAR package.
//
// The CAVIAR package is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the CAVIAR distribution.
//
//========================================================================

#include "caviar2/tools/atom_group.hpp"
#include "caviar2/tools/atom.hpp"
#include "caviar2/tools/atom_list.hpp"

namespace caviar2
{

  namespace tools
  {

    Atom_group::Atom_group(Caviar2 *fptr) : part_of_a_atom_group{false}, upper_level_atom_group{nullptr},
                                            position{Vector3d<double>{0, 0, 0}},
                                            velocity{Vector3d<double>{0, 0, 0}}
    {
    }

    Atom_group::~Atom_group()
    {
    }

    /*
    void Atom_group::verify_settings () {

    }*/

    Atom_group::Atom_group(const Atom_group &a)
    {
    }
    /*
      bool Atom_group::read(caviar2::interpreter::Parser *parser)
      {
        FC_OBJECT_READ_INFO

        while (true)
        {
          FC_IF_RAW_TOKEN_EOF_EOL
          FC_OBJECT_READ_INFO_STR
          FC_IF_GET_REAL3D(position)
          else FC_IF_GET_REAL3D(velocity) else if (string_cmp(t, "add_atom"))
          {
            FIND_OBJECT_BY_NAME(tools, it)
            FC_CHECK_OBJECT_CLASS_NAME(tools, it, atom)
            auto a = *dynamic_cast<tools::Atom *>(object_container->tools[it->second.index]);

            Vector3d<double> pos{0., 0., 0.};
            auto t = parser->get_raw_token();
            std::string ts = t.string_value;
            if (string_cmp(ts, "at_position"))
            {
              GET_OR_CHOOSE_A_REAL_3D_VECTOR(pos, "", "")
            }

            a.position = pos;
            a.upper_level_atom_group = this;
            a.part_of_a_atom_group = true;
            atoms.push_back(a);
            continue;
          }
          else if (string_cmp(t, "clear"))
          {
            atoms.clear();
            continue;
          }
          else FC_ERR_UNDEFINED_VAR(t)
        }

        return true;
      }
    */
    void Atom_group::add_atom(const tools::Atom &a)
    {
      atoms.push_back(a);
    }

    void Atom_group::add_atom(const tools::Atom &a,
                              caviar2::Vector3d<double> p,
                              caviar2::Vector3d<double> v)
    {
      auto at = a;
      at.position = at.position + p;
      at.velocity = at.velocity + v;
      at.upper_level_atom_group = this;
      at.part_of_a_atom_group = true;
      atoms.push_back(at);
    }

    Vector3d<double> Atom_group::pos_tot() const
    {
      if (part_of_a_atom_group)
        return position + upper_level_atom_group->pos_tot();
      else
        return position;
    }

    Vector3d<double> Atom_group::vel_tot() const
    {
      if (part_of_a_atom_group)
        return velocity + upper_level_atom_group->vel_tot();
      else
        return velocity;
    }

  } // tools

}
