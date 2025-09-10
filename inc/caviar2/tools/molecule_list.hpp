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

#include "caviar2/tools/molecule.hpp"

namespace caviar2 {
class Caviar2;
class Atom_data;
namespace tools
{

  /**
   * This class creates list of molecules.
   * list contains references of the molecules.
   */
  class Molecule_list 
  {
  public:
    Molecule_list(class Caviar2* caviar) ;
    Molecule_list(const Molecule_list &);
    Molecule_list();
    ~Molecule_list();

    
    void verify_settings();
    void add_molecule(const tools::Molecule &);
    std::vector<tools::Molecule *> molecules;
        // FC_BASE_OBJECT_COMMON_TOOLS
    void set_caviar(Caviar2 *c) ;

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

} // tools

}

