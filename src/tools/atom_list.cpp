
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

#include "caviar2/tools/atom_list.hpp"
#include "caviar2/tools/atom.hpp"
#include "caviar2/tools/atom_group.hpp"

namespace caviar2
{

  namespace tools
  {

    Atom_list::Atom_list(Caviar2 *fptr)
    {
    }

    Atom_list::~Atom_list()
    {
    }

    void Atom_list::verify_settings()
    {
    }

    Atom_list::Atom_list(const Atom_list &a)
    {
    }
    /*
      bool Atom_list::read(caviar2::interpreter::Parser *parser)
      {
        FC_OBJECT_READ_INFO

        while (true)
        {
          FC_IF_RAW_TOKEN_EOF_EOL
          FC_OBJECT_READ_INFO_STR
          if (string_cmp(t, "add_atom"))
          {
            FIND_OBJECT_BY_NAME(tools, it)
            FC_CHECK_OBJECT_CLASS_NAME(tools, it, atom)
            auto a = dynamic_cast<tools::Atom *>(object_container->tools[it->second.index]);

            atoms.push_back(a);
            continue;
          }
          else if (string_cmp(t, "clear"))
          {
            atoms.clear();
            continue;
          }
          else
            FC_ERR_UNDEFINED_VAR(t)
        }

        return true;
      }
    */
    void Atom_list::add_atom(const tools::Atom &)
    {
    }

  } // tools

}
