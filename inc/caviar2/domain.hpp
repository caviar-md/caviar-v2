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

#include <vector>
#include "caviar2/vector3d.hpp"

namespace caviar2
{

  /**
   * This class is the base class for all the domains.
   * Domain contains the value of the simulation boxes.
   *
   */
  class Domain
  {
  public:
    Domain();
    Domain(class Caviar2* caviar) ;
    virtual ~Domain();

  
  

    virtual void calculate_local_domain();
    virtual void generate();
    virtual void calculate_procs_grid();

    /**
     * It must be called after the changes in the domain
     */
    virtual void update_after_domain_change();

    /**
     * Used in barostat scaling for geometrical forces.
     *
     */
    virtual void scale_position(double scale_ratio, caviar2::Vector3d<int> scale_axis);

    /**
     * calculates process rank from it grid index
     */
    virtual int grid2rank(int x, int y, int z);

    /**
     * with respect to the number of processes, makes a grid with lowest shared area possible.
     */
    virtual void find_best_grid();

    virtual double fix_distance_x(double d);
    virtual double fix_distance_y(double d);
    virtual double fix_distance_z(double d);

    virtual caviar2::Vector3d<double> fix_distance(caviar2::Vector3d<double> v);

    /**
     * Gives the corrected position in periodic boundary condition. Use case: if an atom of a molecule crossed the  global boundary
     * in periodic condition, its actual location is on the other side of the boundary, but since exchanging just one atom of a molecule
     * between domain is MPI expensive, we use fix_position to get the correct location of the atoms in cases it is needed, such as
     * when using the dealii force_field and we need all the particles remain inside mesh.
     */
    virtual caviar2::Vector3d<double> fix_position(caviar2::Vector3d<double> v, caviar2::Vector3d<int> &msd_value, bool &update_verlet_list);

    virtual Vector3d<double> periodic_distance(const Vector3d<double>);

    /**
     * Total volume of (mpi) local domain. In non-mpi simulations, local == global
     */
    virtual double volume_local();

    /**
     * Total volume of global domain
     */
    virtual double volume_global();

    Vector3d<int> boundary_condition;

    int grid_index_x, grid_index_y, grid_index_z; // starts from (0) to (nprocs_i-1) ; i=x,y,z
    int nprocs_x, nprocs_y, nprocs_z;             // it can be at least (1) and at most (nprocs)

    Vector3d<double> lower_global, upper_global;
    Vector3d<double> lower_local, upper_local;

    Vector3d<double> size_local, size_global;

    /**
     * used in MD_MPI case:
     * MPI process rank and number of processes
     */
    int me, nprocs;

    /**
     *  all neighborlist domains around the me=all [1][1][1]. left=all[0][1][1]. up=all[1][2][1].
     *  if one domain exists, in one process case, all[i][j][k]=me for 0<=i,j,k<=2
     *  left&right: x direction, down&up: y direction, bottom&top: z direction. right&up&top are the positive directions
     */
    int all[3][3][3];

    /**
     * defined to have a faster MPI when we have less than 27 processes or when domains can have similar neigbors
     */
    std::vector<int> neighborlist_domains;

    Vector3d<double> half_edge;

    virtual void set_caviar(Caviar2 *c);

    void verify_settings();

  protected:
    Caviar2 *caviar_; // non-owning pointer; or use std::weak_ptr if shared ownership
  };

}