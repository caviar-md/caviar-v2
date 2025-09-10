
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

#include "caviar2/tools/distribution.hpp"
#include "caviar2/atom_data.hpp" // TODO CLEAN THIS PART
#include "caviar2/shape.hpp"
#include "caviar2/tools/grid_1d.hpp"
#include "caviar2/tools/random_1d.hpp"
#include "caviar2/tools/atom.hpp"
#include "caviar2/tools/atom_group.hpp"
#include "caviar2/tools/molecule.hpp"
#include "caviar2/tools/molecule_group.hpp"
#include "caviar2/caviar2.hpp"


namespace caviar2 {

namespace tools
{

  Distribution::Distribution(Caviar2 *fptr) : 
                                             check_radius{false},
                                             atom_data{nullptr}, // TODO CLEAN THIS PART
                                             boundary_shape{nullptr},
                                             atom{nullptr}, atom_group{nullptr},
                                             molecule{nullptr}, molecule_group{nullptr},
                                             grid_1d_x{nullptr}, grid_1d_y{nullptr}, grid_1d_z{nullptr},
                                             random_1d_x{nullptr}, random_1d_y{nullptr}, random_1d_z{nullptr} {
                                                                                             }

                                             Distribution::~Distribution()
  {
  }

  void Distribution::verify_settings()
  {
  }
/*
  bool Distribution::read(caviar2::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;

    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      const auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
      if (string_cmp(t, "boundary_shape"))
      {
        FIND_OBJECT_BY_NAME(shape, it)
        boundary_shape = object_container->shape[it->second.index];
      }
      else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
      {
        FIND_OBJECT_BY_NAME(atom_data, it)
        atom_data = object_container->atom_data[it->second.index];
      }
      else if (string_cmp(t, "atom"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        FC_CHECK_OBJECT_CLASS_NAME(tools, it, atom)
        atom = static_cast<tools::Atom *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "atom_group"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        FC_CHECK_OBJECT_CLASS_NAME(tools, it, atom_group)
        atom_group = static_cast<tools::Atom_group *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "molecule"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        FC_CHECK_OBJECT_CLASS_NAME(tools, it, molecule)
        molecule = static_cast<tools::Molecule *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "molecule_group"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        FC_CHECK_OBJECT_CLASS_NAME(tools, it, molecule_group)
        molecule_group = static_cast<tools::Molecule_group *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "distribute_grid_3d"))
      {
        distribute_grid_3D();
        return in_file;
      }
      else if (string_cmp(t, "distribute_random_3d"))
      {
        int num_of_atoms = 0;
        double radius = 0.0;
        GET_OR_CHOOSE_A_INT(num_of_atoms, "", "")
        GET_OR_CHOOSE_A_REAL(radius, "", "")
        distribute_random_3D(num_of_atoms, radius);
        return in_file;
      }
      else if (string_cmp(t, "grid_1d_x"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        grid_1d_x = static_cast<tools::Grid_1D *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_y"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        grid_1d_y = static_cast<tools::Grid_1D *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "grid_1d_z"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        grid_1d_z = static_cast<tools::Grid_1D *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "random_1d_x"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        random_1d_x = static_cast<tools::Random_1D *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "random_1d_y"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        random_1d_y = static_cast<tools::Random_1D *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "random_1d_z"))
      {
        FIND_OBJECT_BY_NAME(tools, it)
        random_1d_z = static_cast<tools::Random_1D *>(object_container->tools[it->second.index]);
      }
      else if (string_cmp(t, "add_radius"))
      {
        int i = 0;
        double m = 0;
        GET_OR_CHOOSE_A_INT(i, "", "")
        GET_OR_CHOOSE_A_REAL(m, "", "")
        // auto ind = parser->get_positive_int();
        // auto m = parser->get_real();
        if (radius_vector.size() < static_cast<unsigned>(i + 1))
        {
          radius_vector.resize(i + 1);
        }
        radius_vector[i] = m;
        if (m == 0)
          radius_vector[i] = 0.0;
        if (m < 0)
          caviar_->log.warning("you have entered a negative value for radius_vector.");
      }
      else
      {
        FC_ERR_UNDEFINED_VAR(t)
      }
    }
    return in_file;
  }
*/
  bool Distribution::distribute_grid_3D()
  {
    caviar_->log.info("distribute_grid_3D:");


    if (atom == nullptr && molecule == nullptr)
    {
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "(atom==nullptr && molecule==nullptr )");
    }

    if (atom != nullptr)
    {
      if (atom_group == nullptr)
      {
        caviar_->log.error_all(FC_FILE_LINE_FUNC, "(atom_group==nullptr)");
      }
    }

    if (molecule != nullptr)
    {
      if (molecule_group == nullptr)
      {
        caviar_->log.error_all(FC_FILE_LINE_FUNC, "(molecule_group==nullptr)");
      }
    }

    std::vector<Vector3d<double>> p_vector; // total positions of atoms
    if (atom != nullptr)
    {
      p_vector.push_back(atom->pos_tot());
    }

    if (molecule != nullptr)
    {
      for (unsigned int i = 0; i < molecule->atoms.size(); ++i)
      {
        p_vector.push_back(molecule->atoms[i].pos_tot());
      }
    }

    for (unsigned int i = 0; i < grid_1d_x->no_points(); ++i)
    {
      double x = grid_1d_x->give_point(i);
      for (unsigned int j = 0; j < grid_1d_y->no_points(); ++j)
      {
        double y = grid_1d_y->give_point(j);
        for (unsigned int k = 0; k < grid_1d_z->no_points(); ++k)
        {
          double z = grid_1d_z->give_point(k);
          const Vector3d<double> p{x, y, z}, v{0, 0, 0};
          /*
           // simple method - just center is inside condition
           if (boundary_shape->is_inside(p)) {
             if (object_atom_check)
               container_molecule->add_atom (*object_atom, p, v);
             if (object_molecule_check)
               container_molecule->add_molecule (*object_molecule, p, v);
           }
         */

          bool inside_flag = true;
          if (boundary_shape != nullptr)
            for (unsigned int m = 0; m < radius_vector.size(); ++m)
              if (!boundary_shape->is_inside(p + p_vector[m], radius_vector[m]))
              {
                inside_flag = false;
                continue;
              }

          if (inside_flag)
          {
            if (atom != nullptr)
              atom_group->add_atom(*atom, p, v);
            if (molecule != nullptr)
              molecule_group->add_molecule(*molecule, p, v);
          }
        }
      }
    }

    return true;
  }

  bool Distribution::distribute_random_3D(const int num_of_atoms, const double)
  {
    caviar_->log.info("DISTRIBUTE_RANDOM_3D: ");

 

    if (atom == nullptr && molecule == nullptr)
    {
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "(atom==nullptr && molecule==nullptr )");
    }

    if (atom != nullptr)
    {
      if (atom_group == nullptr)
      {
        caviar_->log.error_all(FC_FILE_LINE_FUNC, "(atom_group==nullptr)");
      }
    }

    if (molecule != nullptr)
    {
      if (molecule_group == nullptr)
      {
        caviar_->log.error_all(FC_FILE_LINE_FUNC, "(molecule_group==nullptr)");
      }
    }

    std::vector<Vector3d<double>> p_vector; // total positions of atoms
    if (atom != nullptr)
    {
      p_vector.push_back(atom->pos_tot());
    }

    if (molecule != nullptr)
    {
      for (unsigned int i = 0; i < molecule->atoms.size(); ++i)
      {
        p_vector.push_back(molecule->atoms[i].pos_tot());
      }
    }

    int sum_of_added = 0;
    while (sum_of_added < num_of_atoms)
    {
      double x = random_1d_x->give_value();
      double y = random_1d_y->give_value();
      double z = random_1d_z->give_value();
      const Vector3d<double> p{x, y, z}, v{0, 0, 0};
      /*
       // simple method - just center is inside condition
       if (boundary_shape->is_inside(p)) {
         if (object_atom_check)
           container_molecule->add_atom (*object_atom, p, v);
         if (object_molecule_check)
           container_molecule->add_molecule (*object_molecule, p, v);
       }
     */

      bool inside_flag = true;
      if (boundary_shape != nullptr)
        for (unsigned int m = 0; m < radius_vector.size(); ++m)
          if (!boundary_shape->is_inside(p + p_vector[m], radius_vector[m]))
          {
            inside_flag = false;
            continue;
          }

      if (inside_flag)
      {
        if (atom != nullptr)
        {
          /*
                  bool is_empty = atom_data->empty_of_atoms(p, radius);
                  if (is_empty) {
                    bool added = atom_group->add_atom (*atom, p, v);
                    if (added)
                    sum_of_added++;
                  }
          */
        }

        if (molecule != nullptr)
        {
          caviar_->log.error_all(FC_FILE_LINE_FUNC, "not-implementet yet.");
        }
      }
    }

    return true;
  }

} // tools

}
