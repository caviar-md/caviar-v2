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

#include "caviar2/tools/atom.hpp"

namespace caviar2 {
class Caviar2;
namespace tools
{

  /**
   * This class creates list of atoms.
   * list contains references of the atoms.
   */
  class Atom_list 
  {
  public:
    Atom_list(class Caviar2* caviar) ;
    Atom_list(const Atom_list &);
    Atom_list();
    ~Atom_list();
    
    void verify_settings();
    void add_atom(const tools::Atom &);
    std::vector<tools::Atom *> atoms;
        // FC_BASE_OBJECT_COMMON_TOOLS
    void set_caviar(Caviar2 *c) ;

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // tools

}

