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

#include "caviar2/atom_data.hpp"


#include <string>

namespace caviar2 {

bool Atom_data::add_xyz_data_file(std::string filepath)
{
  // caviar_->log.info("Basic::add_xyz_data_file ");

  // std::string xyz_file_name = "";
  // bool last_frame = false;
  // bool in_file = true;
  // bool read_velocity = false;
  // bool replace_data = false; // used for resuming simulations when there's bonds and angles
  // while (true)
  // {
  //   GET_A_TOKEN_FOR_CREATION
  //   auto t = token.string_value;
  //   if (string_cmp(t, "last_frame"))
  //   {
  //     last_frame = true;
  //   }
  //   else if (string_cmp(t, "replace_data"))
  //   {
  //     replace_data = true;
  //   }
  //   else if (string_cmp(t, "file_name"))
  //   {
  //     const auto token = parser->get_val_token();
  //     const auto file_name = token.string_value;
  //     //      GET_A_STRING(xyz_file_name,"","")
  //     xyz_file_name = file_name;
  //   }
  //   else if (string_cmp(t, "read_velocity"))
  //   {
  //     read_velocity = true;
  //   }
  //   else
  //     FC_ERR_UNDEFINED_VAR(t)
  // }

  // int start_line = 1;
  // std::string xyz_file_name_full = join_path(fptr->input_file_directory, xyz_file_name);
  // if (!file_exists_1(xyz_file_name_full))
  //   caviar_->log.error_all(FC_FILE_LINE_FUNC_PARSE, "file does not exist : " + xyz_file_name_full);

  // if (last_frame)
  // {
  //   int i = 1;
  //   int num_xyz_frames = 0;
  //   //    Parser *pf (xyz_file_name);
  //   caviar2::interpreter::Parser pf(fptr, xyz_file_name_full);
  //   auto t = pf.get_val_token();
  //   while (t.kind != caviar2::interpreter::Kind::eof)
  //   {
  //     if (t.kind == caviar2::interpreter::Kind::identifier)
  //     {
  //       auto ts = t.string_value;
  //       if (string_cmp(ts, "atom") || string_cmp(ts, "atoms") || string_cmp(ts, "ATOM") || string_cmp(ts, "ATOMS"))
  //       {
  //         ++num_xyz_frames;
  //         if (i == 1)
  //           caviar_->log.error_all(FC_FILE_LINE_FUNC_PARSE, "Unknown xyz format");
  //         start_line = i - 1;
  //       }
  //     }
  //     pf.end_of_line();
  //     t = pf.get_val_token();
  //     ++i;
  //   }
  //   std::string st1 = "Num. of total xyz frames found : " + std::to_string(num_xyz_frames);
  //   std::string st2 = "Last xyz frame found at line number : " + std::to_string(i);
  //   caviar_->log.info(st1);
  //   caviar_->log.info(st2);
  // }

  // caviar2::interpreter::Parser pf(fptr, xyz_file_name_full);
  // for (int i = 1; i < start_line; ++i)
  // {
  //   pf.end_of_line();
  // }

  // /*
  //   int num_atoms = pf.get_literal_int();
  //   std::cout << "num_atoms : " << num_atoms << std::endl;
  //       pf.end_of_line();
  //   auto st = pf.get_val_token();
  //   std::cout << "st : " << st.string_value << std::endl;
  //       pf.end_of_line();
  // */

  // int num_atoms = pf.get_literal_int();

  // pf.end_of_line();

  // std::string st = " Number of atoms found at the inputted xyz file : " + std::to_string(num_atoms);

  // if (num_atoms == 0)
  //   caviar_->log.error_all(FC_FILE_LINE_FUNC, "no atom found in the xyz file. Something is wrong.");

  // caviar_->log.info(st);
  // pf.get_val_token();
  // pf.end_of_line();

  // if (replace_data)
  // {
  //   if (num_atoms != (int)atom_struct_owned.position.size())
  //     caviar_->log.error_all(FC_FILE_LINE_FUNC, "XYZ file is not compatible: Different number of existing atoms in atomdata and xyz file atoms.");
  // }

  // for (int i = 0; i < num_atoms; ++i)
  // {
  //   auto type = pf.get_literal_int();

  //   Vector3d<double> pos, vel{0.0, 0.0, 0.0};
  //   ;

  //   pos.x = pf.get_literal_real();
  //   pos.y = pf.get_literal_real();
  //   pos.z = pf.get_literal_real();

  //   if (read_velocity)
  //   {
  //     vel.x = pf.get_literal_real();
  //     vel.y = pf.get_literal_real();
  //     vel.z = pf.get_literal_real();
  //   }

  //   pf.end_of_line();

  //   if (replace_data)
  //   {
  //     if (type != (int)atom_struct_owned.type[i])
  //       caviar_->log.error_all(FC_FILE_LINE_FUNC, "XYZ file is not compatible: Different atom type order exists.");
  //     atom_struct_owned.position[i] = pos;
  //     atom_struct_owned.velocity[i] = vel;
  //   }
  //   else
  //   {
  //     auto id = get_num_of_atoms_global();
  //     add_atom(id, type, pos, vel);
  //   }
  // }

  //return in_file;
  return true;
}

}
