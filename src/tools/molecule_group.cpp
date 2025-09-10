
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

#include "caviar2/tools/molecule_group.hpp"
#include "caviar2/tools/molecule.hpp"
#include "caviar2/tools/molecule_list.hpp"

namespace caviar2
{

  namespace tools
  {

    Molecule_group::Molecule_group(Caviar2 *fptr) : part_of_a_molecule_group{false}, upper_level_molecule_group{nullptr},
                                                    position{Vector3d<double>{0, 0, 0}},
                                                    velocity{Vector3d<double>{0, 0, 0}}
    {
    }

    Molecule_group::~Molecule_group()
    {
    }

    void Molecule_group::verify_settings()
    {
    }

    Molecule_group::Molecule_group(const Molecule_group &a) 
    {
    }
    /*
      bool Molecule_group::read(caviar2::interpreter::Parser *parser)
      {
        FC_OBJECT_READ_INFO

        while (true)
        {
          FC_IF_RAW_TOKEN_EOF_EOL
          FC_OBJECT_READ_INFO_STR
          FC_IF_GET_REAL3D(position)
          else FC_IF_GET_REAL3D(velocity) else if (string_cmp(t, "add_molecule"))
          {
            FIND_OBJECT_BY_NAME(tools, it)
            FC_CHECK_OBJECT_CLASS_NAME(tools, it, molecule)
            auto m = *dynamic_cast<tools::Molecule *>(object_container->tools[it->second.index]);

            Vector3d<double> pos{0., 0., 0.};
            auto token = parser->get_raw_token();
            std::string t = token.string_value;
            if (string_cmp(t, "at_position"))
            {
              GET_OR_CHOOSE_A_REAL_3D_VECTOR(pos, "", "")
            }

            m.position = pos;
            m.upper_level_molecule_group = this;
            m.part_of_a_molecule_group = true;
            molecules.push_back(m);
            continue;
          }
          else if (string_cmp(t, "clear"))
          {
            molecules.clear();
            continue;
          }
          else FC_ERR_UNDEFINED_VAR(t)
        }

        return true;
      }
    */
    void Molecule_group::add_molecule(const tools::Molecule &m)
    {
      molecules.push_back(m);
    }

    void Molecule_group::add_molecule(const tools::Molecule &m,
                                      caviar2::Vector3d<double> p,
                                      caviar2::Vector3d<double> v)
    {
      auto mt = m;
      mt.position = mt.position + p;
      mt.velocity = mt.velocity + v;
      mt.upper_level_molecule_group = this;
      mt.part_of_a_molecule_group = true;
      // std::cout << mt.position << std::endl;
      molecules.push_back(mt);
    }

    Vector3d<double> Molecule_group::pos_tot() const
    {
      if (part_of_a_molecule_group)
        return position + upper_level_molecule_group->pos_tot();
      else
        return position;
    }

    Vector3d<double> Molecule_group::vel_tot() const
    {
      if (part_of_a_molecule_group)
        return velocity + upper_level_molecule_group->vel_tot();
      else
        return velocity;
    }

  } // tools

}
