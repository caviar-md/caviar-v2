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

#include "caviar2/vector3d.hpp"

namespace caviar2
{

    class Caviar2;
    class Atom_data;
    class Domain;
    class Neighborlist;

    /**
     * @class Force_field
     * @brief Base class for all force fields in the Caviar2 library.
     *
     * This is an abstract base class providing the interface for derived force field classes.
     * Derived classes must implement all pure virtual methods.
     *
     * It manages key components like atoms, domain, neighbor list, and handles energy and force calculations.
     */
    class Force_field
    {
    public:
        /**
         * @brief Constructor.
         * @param caviar Pointer to the main Caviar object.
         */
        Force_field(class Caviar2 *caviar);

        Force_field();

        virtual void set_caviar(Caviar2 *c);

        /**
         * @brief Virtual destructor.
         */
        virtual ~Force_field() = default;

        /**
         * @brief Calculates accelerations/forces on atoms.
         *
         * Must be implemented by derived classes.
         */
        virtual void calculate_acceleration() = 0;

        /**
         * @brief Calculates the total energy of the system under this force field.
         *
         * @return Total energy as a double.
         */
        virtual double energy();

        /**
         * @brief Calculates potential energy for a given vector.
         *
         * @param vec Input vector (e.g. position or displacement).
         * @return Potential energy as a double.
         */
        virtual double potential(const Vector3d<double> &vec);

        /**
         * @brief Calculates potential energy for an atom specified by its index.
         *
         * @param index Index of the atom.
         * @return Potential energy as a double.
         */
        virtual double potential(int index);

        /**
         * @brief Calculates the force field vector for a given vector.
         *
         * @param vec Input vector.
         * @return Force field vector.
         */
        virtual Vector3d<double> field(const Vector3d<double> &vec);

        /**
         * @brief Calculates the force field vector for an atom specified by its index.
         *
         * @param index Index of the atom.
         * @return Force field vector.
         */
        virtual Vector3d<double> field(int index);

        /**
         * @brief Scales atomic positions for barostat or other geometric transformations.
         *
         * @param scaleRatio The scaling ratio applied to positions.
         * @param scaleAxis The axis along which scaling is applied.
         */
        virtual void scale_position(double scaleRatio, Vector3d<int> scaleAxis);

        /** @brief Cutoff distance for the force calculations. */
        double cutoff = 0.0;

        /** @brief Pointer to atom data used in the force field. */
        Atom_data *atom_data = nullptr;

        /** @brief Pointer to the simulation domain. */
        Domain *domain = nullptr;

        /** @brief Pointer to the neighbor list for efficient neighbor searches. */
        Neighborlist *neighborlist = nullptr;

        /** @brief MPI rank for parallel computations. */
        int my_mpi_rank = -1;

        virtual void verify_settings();


    protected:
        Caviar2 *caviar_ = nullptr; // non-owning pointer; or use std::weak_ptr if shared ownership
    };

} // namespace caviar2