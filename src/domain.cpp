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
#include "caviar2/caviar2.hpp"

#include "caviar2/domain.hpp"
#include "caviar2/communicator.hpp"

#include <array>
#ifdef CAVIAR_WITH_MPI
#include <mpi.h>
#endif

namespace caviar2
{
  Domain::Domain() {}

  Domain::Domain(Caviar2 *fptr) : caviar_{fptr},
                                  boundary_condition{Vector3d<int>{0, 0, 0}},
                                  grid_index_x{0}, grid_index_y{0}, grid_index_z{0},
                                  nprocs_x{1}, nprocs_y{1}, nprocs_z{1}, me{0}, nprocs{1}
  {

#if defined(CAVIAR_WITH_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
  }

  Domain::~Domain()
  {
  }

  void Domain::set_caviar(caviar2::Caviar2 *c)
  {
    caviar_ = c;
  }

  void Domain::verify_settings()
  {
  }

  int Domain::grid2rank(int x, int y, int z)
  {
    return x + y * nprocs_x + z * nprocs_x * nprocs_y;
  } // calculates process rank from it grid index

  void Domain::calculate_procs_grid()
  {
#ifdef CAVIAR_WITH_MPI
    neighborlist_domains.push_back(me); // don't know if it's neccesary
    find_best_grid();
    if (me == 0)
      std::cout << "[INF] Number of processes:\n\tnpx:" << nprocs_x << " npy:"
                << nprocs_y << " npz:" << nprocs_z << std::endl;

    grid_index_x = me % nprocs_x;
    grid_index_y = me / nprocs_x % nprocs_y;
    grid_index_z = me / nprocs_x / nprocs_y;

    for (auto i = -1; i < 2; ++i)
      for (auto j = -1; j < 2; ++j)
        for (auto k = -1; k < 2; ++k)
        {
          int ii = 0, jj = 0, kk = 0; // Zero value is just to stop '-Wmaybe-uninitialized'
          if (i == 0)
            ii = grid_index_x;
          if (i == -1)
            ii = grid_index_x == 0 ? nprocs_x - 1 : (grid_index_x - 1) % nprocs_x;
          if (i == +1)
            ii = grid_index_x == nprocs_x - 1 ? 0 : (grid_index_x + 1) % nprocs_x;
          if (j == 0)
            jj = grid_index_y;
          if (j == -1)
            jj = grid_index_y == 0 ? nprocs_y - 1 : (grid_index_y - 1) % nprocs_y;
          if (j == +1)
            jj = grid_index_y == nprocs_y - 1 ? 0 : (grid_index_y + 1) % nprocs_y;
          if (k == 0)
            kk = grid_index_z;
          if (k == -1)
            kk = grid_index_z == 0 ? nprocs_z - 1 : (grid_index_z - 1) % nprocs_z;
          if (k == +1)
            kk = grid_index_z == nprocs_z - 1 ? 0 : (grid_index_z + 1) % nprocs_z;
          all[i + 1][j + 1][k + 1] = grid2rank(ii, jj, kk);

          bool domain_found = false;
          for (auto l : neighborlist_domains)
            if (all[i + 1][j + 1][k + 1] == l)
              domain_found = true;
          if (!domain_found)
            neighborlist_domains.push_back(all[i + 1][j + 1][k + 1]);
        }
#else
    neighborlist_domains.push_back(0); // don't know if it's neccesary
#endif
  }

#if defined(CAVIAR_WITH_MPI)
  static std::vector<std::array<int, 3>> possible_grids(int nprocs)
  {
    std::vector<std::array<int, 3>> grids;
    for (auto nprocs_x = 1; nprocs_x <= nprocs; ++nprocs_x)
    {
      if (nprocs % nprocs_x)
        continue;
      auto nprocs_yz = nprocs / nprocs_x;
      for (auto nprocs_y = 1; nprocs_y <= nprocs_yz; ++nprocs_y)
      {
        if (nprocs_yz % nprocs_y)
          continue;
        // std::array<int, 3> grid = {nprocs_x, nprocs_y, nprocs_yz / nprocs_y};
        // grids.push_back(grid);
        grids.push_back({nprocs_x, nprocs_y, nprocs_yz / nprocs_y});
      }
    }
    return grids;
  }
#endif

  void Domain::find_best_grid()
  {
#if defined(CAVIAR_WITH_MPI)
    auto box_length_x = upper_global.x - lower_global.x;
    auto box_length_y = upper_global.y - lower_global.y;
    auto box_length_z = upper_global.z - lower_global.z;

    auto area_xy = box_length_x * box_length_y, area_xz = box_length_x * box_length_z, area_yz = box_length_y * box_length_z;

    auto grids = possible_grids(nprocs);
    auto min_area = 1.1 * (area_xy + area_xz + area_yz);
    auto min_area_index = 0;
    for (unsigned int i = 0; i < grids.size(); ++i)
    {
      auto area = area_xy / grids[i][0] / grids[i][1] + area_xz / grids[i][0] / grids[i][2] + area_yz / grids[i][1] / grids[i][2];
      if (area < min_area)
      {
        min_area = area;
        min_area_index = i;
      }
    }
    nprocs_x = grids[min_area_index][0];
    nprocs_y = grids[min_area_index][1];
    nprocs_z = grids[min_area_index][2];
#endif
  }

  // note that this force cannot be used when one is using neighborlist loop.
  // currently its main usage is when one is using Shake like algorithms.
  // or Spring_bond or Spring_angle force_fields.
  Vector3d<double> Domain::periodic_distance(const Vector3d<double> v)
  {
    caviar2::Vector3d<double> vf = v;
    static caviar2::Vector3d<double> domain_dh = {0.5 * (upper_global.x - lower_global.x),
                                                  0.5 * (upper_global.y - lower_global.y),
                                                  0.5 * (upper_global.z - lower_global.z)};
    if (boundary_condition.x == 1)
    {
      while (vf.x > +domain_dh.x)
      {
        vf.x -= domain_dh.x * 2.0;
      }
      while (vf.x < -domain_dh.x)
      {
        vf.x += domain_dh.x * 2.0;
      }
    }

    if (boundary_condition.y == 1)
    {
      while (vf.y > +domain_dh.y)
      {
        vf.y -= domain_dh.y * 2.0;
      }
      while (vf.y < -domain_dh.y)
      {
        vf.y += domain_dh.y * 2.0;
      }
    }

    if (boundary_condition.z == 1)
    {
      while (vf.z > +domain_dh.z)
      {
        vf.z -= domain_dh.z * 2.0;
      }
      while (vf.z < -domain_dh.z)
      {
        vf.z += domain_dh.z * 2.0;
      }
    }

    // std::cout << domain_dh << std::endl;
    // if( !( v.x == vf.x && v.y == vf.y && v.z == vf.z) )
    // std::cout << "x: " << v << " vs " << vf << std::endl;
    return vf;
  }

  /*
  bool Domain::read(caviar2::interpreter::Parser *parser)
  {
    FC_OBJECT_READ_INFO
    bool in_file = true;
    while (true)
    {
      GET_A_TOKEN_FOR_CREATION
      auto t = token.string_value;
      FC_OBJECT_READ_INFO_STR
      if (string_cmp(t, "lower_global.x") || string_cmp(t, "xmin"))
      {
        GET_OR_CHOOSE_A_REAL(lower_global.x, "", "")
      }
      else if (string_cmp(t, "upper_global.x") || string_cmp(t, "xmax"))
      {
        GET_OR_CHOOSE_A_REAL(upper_global.x, "", "")
      }
      else if (string_cmp(t, "lower_global.y") || string_cmp(t, "ymin"))
      {
        GET_OR_CHOOSE_A_REAL(lower_global.y, "", "")
      }
      else if (string_cmp(t, "upper_global.y") || string_cmp(t, "ymax"))
      {
        GET_OR_CHOOSE_A_REAL(upper_global.y, "", "")
      }
      else if (string_cmp(t, "lower_global.z") || string_cmp(t, "zmin"))
      {
        GET_OR_CHOOSE_A_REAL(lower_global.z, "", "")
      }
      else if (string_cmp(t, "upper_global.z") || string_cmp(t, "zmax"))
      {
        GET_OR_CHOOSE_A_REAL(upper_global.z, "", "")
      }
      else if (string_cmp(t, "boundary_condition") || string_cmp(t, "bc"))
      {
        GET_OR_CHOOSE_A_INT_3D_VECTOR(boundary_condition, "", "")
      }
      else if (string_cmp(t, "me"))
      {
        std::cout << "ME: " << me << std::endl;
      }
      else if (string_cmp(t, "info"))
      {
        std::cout << "MPI Rank: " << me
                  << " local.x [" << lower_local.x << " , " << upper_local.x << "]"
                  << " local.y [" << lower_local.y << " , " << upper_local.y << "]"
                  << " local.z [" << lower_local.z << " , " << upper_local.z << "]"
                  << " global.x [" << lower_global.x << " , " << upper_global.x << "]"
                  << " global.y [" << lower_global.y << " , " << upper_global.y << "]"
                  << " global.z [" << lower_global.z << " , " << upper_global.z << "]"
                  << std::endl;
      }
      else if (string_cmp(t, "generate"))
      {
        generate();
      }
      else
        caviar_->log.error_all(FC_FILE_LINE_FUNC_PARSE, "Unknown variable or command");
    }
    return in_file;
  }
  */

  void Domain::set_lower_global(const Vector3d<double> &v)
  {
    lower_global = v;
  }
  void Domain::set_upper_global(const Vector3d<double> &v)
  {
    upper_global = v;
  }
  void Domain::set_boundary_condition(const Vector3d<int> &v)
  {
    boundary_condition = v;
  }
  void Domain::print_info()
  {
    std::cout << "MPI Rank: " << me
              << " local.x [" << lower_local.x << " , " << upper_local.x << "]"
              << " local.y [" << lower_local.y << " , " << upper_local.y << "]"
              << " local.z [" << lower_local.z << " , " << upper_local.z << "]"
              << " global.x [" << lower_global.x << " , " << upper_global.x << "]"
              << " global.y [" << lower_global.y << " , " << upper_global.y << "]"
              << " global.z [" << lower_global.z << " , " << upper_global.z << "]"
              << std::endl;
  }
  void Domain::generate()
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
    calculate_local_domain();
#else
    calculate_procs_grid();
    calculate_local_domain();
    std::string s = "boundary condition: ";
    s += std::to_string(boundary_condition.x) + " " + std::to_string(boundary_condition.y) + " " + std::to_string(boundary_condition.z);
    caviar_->log.info(s);
#endif
    update_after_domain_change();
  }

  void Domain::calculate_local_domain()
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
    lower_local.x = lower_global.x;
    lower_local.y = lower_global.y;
    lower_local.z = lower_global.z;

    upper_local.x = upper_global.x;
    upper_local.y = upper_global.y;
    upper_local.z = upper_global.z;
#else
    lower_local.x = lower_global.x + (upper_global.x - lower_global.x) * grid_index_x / nprocs_x;
    lower_local.y = lower_global.y + (upper_global.y - lower_global.y) * grid_index_y / nprocs_y;
    lower_local.z = lower_global.z + (upper_global.z - lower_global.z) * grid_index_z / nprocs_z;

    upper_local.x = lower_global.x + (upper_global.x - lower_global.x) * (grid_index_x + 1) / nprocs_x;
    upper_local.y = lower_global.y + (upper_global.y - lower_global.y) * (grid_index_y + 1) / nprocs_y;
    upper_local.z = lower_global.z + (upper_global.z - lower_global.z) * (grid_index_z + 1) / nprocs_z;
#endif
  }

  double Domain::fix_distance_x(double d)
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
    return d;
#endif
    if (boundary_condition.x == 1)
    {
      if (d > +half_edge.x)
      {
        d -= half_edge.x * 2.0;
      }
      if (d < -half_edge.x)
      {
        d += half_edge.x * 2.0;
      }
    }
    return d;
  }

  double Domain::fix_distance_y(double d)
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
    return d;
#endif
    if (boundary_condition.y == 1)
    {
      if (d > +half_edge.y)
      {
        d -= half_edge.y * 2.0;
      }
      if (d < -half_edge.y)
      {
        d += half_edge.y * 2.0;
      }
    }
    return d;
  }

  double Domain::fix_distance_z(double d)
  {
#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

#elif defined(CAVIAR_WITH_MPI)
    return d;
#endif
    if (boundary_condition.z == 1)
    {
      if (d > +half_edge.z)
      {
        d -= half_edge.z * 2.0;
      }
      if (d < -half_edge.z)
      {
        d += half_edge.z * 2.0;
      }
    }
    return d;
  }

  caviar2::Vector3d<double> Domain::fix_distance(caviar2::Vector3d<double> v)
  {
    return caviar2::Vector3d<double>{fix_distance_x(v.x),
                                     fix_distance_y(v.y),
                                     fix_distance_z(v.z)};
  }

  void Domain::scale_position(double xi, caviar2::Vector3d<int> scale_axis)
  {
    bool x_axis = (scale_axis.x == 1 ? true : false);
    bool y_axis = (scale_axis.x == 1 ? true : false);
    bool z_axis = (scale_axis.x == 1 ? true : false);

    if (x_axis)
    {
      lower_global.x *= xi;
      upper_global.x *= xi;
    }

    if (y_axis)
    {
      lower_global.y *= xi;
      upper_global.y *= xi;
    }

    if (z_axis)
    {
      lower_global.z *= xi;
      upper_global.z *= xi;
    }

    calculate_local_domain();

    update_after_domain_change();

    bool debug = false;
    if (debug)
      std::cout << "doma: [" << lower_global.x << " , " << upper_global.x << "] , [" << lower_global.y << " , " << upper_global.y << "] , [" << lower_global.z << " , " << upper_global.z << "]" << std::endl;
  }

  double Domain::volume_local()
  {
    return size_local.x * size_local.y * size_local.z;
  }

  double Domain::volume_global()
  {
    return size_global.x * size_global.y * size_global.z;
  }

  caviar2::Vector3d<double> Domain::fix_position(caviar2::Vector3d<double> p, caviar2::Vector3d<int> &msd, bool &update_verlet_list)
  {

    if (boundary_condition.x == 1)
    {
      if (p.x < lower_global.x) // while or if
      {
        update_verlet_list = true;
        p.x += size_global.x;
        msd.x -= 1;
      }
      else if (p.x > upper_global.x)
      {
        update_verlet_list = true;
        p.x -= size_global.x;
        msd.x += 1;
      }
    }
    if (boundary_condition.y == 1)
    {
      if (p.y < lower_global.y)
      {
        update_verlet_list = true;
        p.y += size_global.y;
        msd.y -= 1;
      }
      else if (p.y > upper_global.y)
      {
        update_verlet_list = true;
        p.y -= size_global.y;
        msd.y += 1;
      }
    }
    if (boundary_condition.z == 1)
    {
      if (p.z < lower_global.z)
      {
        update_verlet_list = true;
        p.z += size_global.z;
        msd.z -= 1;
      }
      else if (p.z > upper_global.z)
      {
        p.z -= size_global.z;
        update_verlet_list = true;
        msd.z += 1;
      }
    }

    return p;
  }

  void Domain::update_after_domain_change()
  {
    size_local = upper_local - lower_local;
    size_global = upper_global - lower_global;

    // it is defined for one-domain case.
    half_edge = 0.5 * (upper_global - lower_global); // XXX
  }

}
