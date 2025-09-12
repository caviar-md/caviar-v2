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
#include "caviar2/md_simulator.hpp"
#include "caviar2/atom_data.hpp"
#include "caviar2/neighborlist.hpp"
#include "caviar2/force_field.hpp"
#include "caviar2/constraint.hpp"
#include "caviar2/writer.hpp"
#include "caviar2/tools/time_function.hpp"
#include "caviar2/tools/time_function_3d.hpp"
#include "caviar2/communicator.hpp"
#include "caviar2/time_utility.hpp"

// #ifdef CAVIAR_WITH_MPI
// #include <mpi.h>
// #endif
namespace caviar2
{

  Md_simulator::Md_simulator() {}

  Md_simulator::Md_simulator(Caviar2 *fptr) : caviar_{fptr},
                                              atom_data{nullptr},
                                              initialized{false}
  {

    current_step = 0;
    use_time = false;
    use_step = false;

    dt = -1.0;

    auto &L = langevin_param;
    L.rnd_tools_x.seed(1);
    L.rnd_tools_y.seed(2);
    L.rnd_tools_z.seed(3);
    L.kb = -1;
    L.temperature = -1;
    L.kbt = -1;
    L.friction = 0.0;
    set_caviar(fptr);
  }

  Md_simulator::~Md_simulator() {}

  void Md_simulator::set_caviar(caviar2::Caviar2 *c)
  {
    caviar_ = c;
    atom_data = &(c->atom_data);
    // force_field = &c->forces;
    // neighborlist = &c->neighborlist;
    // constraint = &c->constraints;
    // writer = &(c->writers);
  }

  void Md_simulator::set_initial_time(double t)
  {
    initial_time = t;
    use_time = true;
  }

  void Md_simulator::set_final_time(double t)
  {
    final_time = t;
    use_time = true;
  }
  void Md_simulator::set_initial_step(size_t i)
  {
    initial_step = i;
    use_step = true;
  }

  void Md_simulator::set_final_step(size_t i)
  {
    final_step = i;
    use_step = true;
  }

  void Md_simulator::set_final_add(size_t i)
  {

    final_step += i;
    initial_step = current_step;
  }
  void Md_simulator::langevin_temperature(double langevin_temperature)
  {
    langevin_param.temperature = langevin_temperature;
    if (langevin_param.temperature < 0.0)
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "Temperature have to non-negative.");
  }
  void Md_simulator::langevin_friction(double langevin_friction)
  {
    langevin_param.friction = langevin_friction;
    if (langevin_param.friction < 0.0)
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "friction have to non-negative.");
  }
  void Md_simulator::langevin_kb(double langevin_kb)
  {
    langevin_param.kb = langevin_kb;
    if (langevin_param.kb < 0.0)
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "kb have to non-negative.");
  }
  void Md_simulator::langevin_kbt(double langevin_kbt)
  {
    langevin_param.kbt = langevin_kbt;
    if (langevin_param.kbt < 0.0)
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "kbt have to non-negative.");
  }

  void Md_simulator::set_integrator_type(Integrator_t i_type)
  {
    integrator_type = i_type;

    if (integrator_type == Integrator_t::Velocity_verlet)
      atom_data->set_record_owned_acceleration_old(true);

    if (integrator_type == Integrator_t::Unknown)
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "invalid or non-implemented integrator");
  }

  void Md_simulator::set_dt(double dt_)
  {
    dt = dt_;
  }

  /*
    bool Md_simulator::read(caviar2::interpreter::Parser *parser)
    {
      FC_OBJECT_READ_INFO
      bool in_file = true;
      while (true)
      {
        GET_A_TOKEN_FOR_CREATION
        auto t = token.string_value;
        FC_OBJECT_READ_INFO_STR
        if (string_cmp(t, "dt"))
        {
          GET_OR_CHOOSE_A_REAL(dt, "", "")
        }
        else if (string_cmp(t, "initial_time"))
        {
          GET_OR_CHOOSE_A_REAL(initial_time, "", "")
          use_time = true;
        }
        else if (string_cmp(t, "final_time"))
        {
          GET_OR_CHOOSE_A_REAL(final_time, "", "")
          use_time = true;
        }
        else if (string_cmp(t, "initial_step"))
        {
          GET_OR_CHOOSE_A_INT(initial_step, "", "")
          use_step = true;
        }
        else if (string_cmp(t, "final_step"))
        {
          GET_OR_CHOOSE_A_INT(final_step, "", "")
          use_step = true;
        }
        else if (string_cmp(t, "final_step_add"))
        {
          int add_step = 0;
          GET_OR_CHOOSE_A_INT(add_step, "", "")
          final_step += add_step;
          initial_step = current_step;
        }
        else if (string_cmp(t, "run"))
        {
          run();
        }
        else if (string_cmp(t, "add_constraint") || string_cmp(t, "constraint"))
        {
          FIND_OBJECT_BY_NAME(constraint, it)
          constraint.push_back(object_container->constraint[it->second.index]);
        }
        else if (string_cmp(t, "remove_constraint"))
        {
  #define FC_REMOVE_OBJECT(name)                               \
    FIND_OBJECT_BY_NAME(name, it)                              \
    bool found = false;                                        \
    for (unsigned int i = 0; i < name.size(); ++i)             \
    {                                                          \
      if (name[i] == object_container->name[it->second.index]) \
      {                                                        \
        name.erase(name.begin() + i);                          \
        found = true;                                          \
        break;                                                 \
      }                                                        \
    }                                                          \
    if (!found)                                                \
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "Unknown " #name " to remove.");
          FC_REMOVE_OBJECT(constraint)
        }
        else if (string_cmp(t, "add_writer") || string_cmp(t, "writer"))
        {
          FIND_OBJECT_BY_NAME(writer, it)
          writer.push_back(object_container->writer[it->second.index]);
        }
        else if (string_cmp(t, "add_neighborlist") || string_cmp(t, "neighborlist"))
        {
          FIND_OBJECT_BY_NAME(neighborlist, it)
          neighborlist.push_back(object_container->neighborlist[it->second.index]);
        }
        else if (string_cmp(t, "add_force_field") || string_cmp(t, "force_field"))
        {
          FIND_OBJECT_BY_NAME(force_field, it)
          force_field.push_back(object_container->force_field[it->second.index]);
        }
        else if (string_cmp(t, "remove_writer"))
        {
          FC_REMOVE_OBJECT(writer)
        }
        else if (string_cmp(t, "remove_neighborlist"))
        {
          FC_REMOVE_OBJECT(neighborlist)
        }
        else if (string_cmp(t, "remove_force_field"))
        {
          FC_REMOVE_OBJECT(force_field)
        }
        else if (string_cmp(t, "set_atom_data") || string_cmp(t, "atom_data"))
        {
          FIND_OBJECT_BY_NAME(atom_data, it)
          atom_data = object_container->atom_data[it->second.index];
        }
        else if (string_cmp(t, "add_time_function") || string_cmp(t, "time_function"))
        {
          FIND_OBJECT_BY_NAME(tools, it)
          FC_CHECK_OBJECT_CLASS_NAME(tools, it, time_function)
          tools::Time_function *a = dynamic_cast<tools::Time_function *>(object_container->tools[it->second.index]);
          time_function.push_back(a);
        }
        else if (string_cmp(t, "add_time_function_3d") || string_cmp(t, "time_function_3d"))
        {
          FIND_OBJECT_BY_NAME(tools, it)
          FC_CHECK_OBJECT_CLASS_NAME(tools, it, time_function_3d)
          tools::Time_function_3d *a = dynamic_cast<tools::Time_function_3d *>(object_container->tools[it->second.index]);
          time_function_3d.push_back(a);
        }
        else if (string_cmp(t, "step"))
        {
          // int i=0;
          // GET_OR_CHOOSE_A_INT(i,"","")
          ++current_step;
          step();
          continue;
        }
        else if (string_cmp(t, "integrator_type"))
        {

          auto t2 = parser->get_val_token();
          auto ts = t2.string_value;
          if (string_cmp(ts, "leap_frog"))
            integrator_type = Integrator_t::Leap_frog;
          else if (string_cmp(ts, "velocity_verlet"))
          {
            if (atom_data == nullptr)
              caviar_->log.error_all(FC_FILE_LINE_FUNC, "set atom_data before this integrator");
            atom_data->set_record_owned_acceleration_old(true);
            integrator_type = Integrator_t::Velocity_verlet;
          }
          else if (string_cmp(ts, "velocity_verlet_langevin"))
            integrator_type = Integrator_t::Velocity_verlet_langevin;
          else
            caviar_->log.error_all(FC_FILE_LINE_FUNC, "invalid or non-implemented integrator");
        }
        else if (string_cmp(t, "temperature"))
        {
          GET_OR_CHOOSE_A_REAL(langevin_param.temperature, "", "")
          if (langevin_param.temperature < 0.0)
            caviar_->log.error_all(FC_FILE_LINE_FUNC, "Temperature have to non-negative.");
        }
        else if (string_cmp(t, "friction"))
        {
          GET_OR_CHOOSE_A_REAL(langevin_param.friction, "", "")
          if (langevin_param.friction < 0.0)
            caviar_->log.error_all(FC_FILE_LINE_FUNC, "friction have to non-negative.");
        }
        else if (string_cmp(t, "kb"))
        {
          GET_OR_CHOOSE_A_REAL(langevin_param.kb, "", "")
          if (langevin_param.kb < 0.0)
            caviar_->log.error_all(FC_FILE_LINE_FUNC, "kb have to non-negative.");
        }
        else if (string_cmp(t, "kbt"))
        {
          GET_OR_CHOOSE_A_REAL(langevin_param.kbt, "", "")
          if (langevin_param.kb < 0.0)
            caviar_->log.error_all(FC_FILE_LINE_FUNC, "kbt have to non-negative.");
        }
        else if (string_cmp(t, "dt"))
        {
          GET_OR_CHOOSE_A_REAL(dt, "", "")
        }
        else if (read_base_class_commands(parser))
        {
        }
        else
          caviar_->log.error_all(FC_FILE_LINE_FUNC, "Unknown variable or command");
      }
      return in_file;
    }
    */
  bool Md_simulator::run()
  {
    caviar_->log.info("MD started.");
#ifdef CAVIAR_WITH_OPENMP
    caviar_->log.info("CAVIAR is started with OpenMP. Max number of threads is " + std::to_string(omp_get_max_threads()));
#endif
    verify_settings();

    double wall0 = get_wall_time();
    // double cpu0  = get_cpu_time();

    // t_start = clock();

    for (current_step = initial_step; current_step < final_step; ++current_step)
      step();

    double wall1 = get_wall_time();
    // double cpu1  = get_cpu_time();

    cleanup();

    caviar_->log.info("MD finished.");

    // t_end = clock();
    // double simulation_time =  ( (float)t_end - (float)t_start ) / CLOCKS_PER_SEC;

    // #if defined (CAVIAR_WITH_MPI)
    //   std::string s = "MPI process number" + std::to_string(comm->me) + " : ";
    // #else
    //   std::string s = "";
    // #endif

    std::string s = "md_simulation run call:";
    s += "Wall Time = " + std::to_string(wall1 - wall0) + " sec. , ";
    // s += "CPU Time  = " + std::to_string(cpu1  - cpu0) + " sec. ";

    caviar_->log.info(s);

    return true; // WARNING
  }

  void Md_simulator::verify_settings()
  {
    if (use_time && use_step)
      caviar_->log.error_all("simulator md: verify_settings: cannot use both 'time' and 'step'.");
    if (!(use_time || use_step))
      caviar_->log.error_all("simulator md: verify_settings: please assign 'time' or 'step' variables.");
  }

  void Md_simulator::initialize()
  {

    atom_data->set_atoms_mpi_rank();

    my_mpi_rank = atom_data->get_mpi_rank();

    if (dt < 0.0)
      caviar_->log.error_all(FC_FILE_LINE_FUNC, "expected 'dt' input");

    switch (integrator_type)
    {

      // case Integrator_t::Verlet :
      // caviar_->log.error_all(FC_FILE_LINE_FUNC,"not implemented");
      // integrate_verlet();

      break;

    default:
    case Integrator_t::Velocity_verlet:
      // integrate_velocity_verlet();
      break;

    case Integrator_t::Leap_frog:
      // integrate_leap_frog();
      break;

    case Integrator_t::Velocity_verlet_langevin:

      auto &L = langevin_param;

      L.a = (2.0 - L.friction * dt) / (2.0 + L.friction * dt);
      if (L.kbt > 0)
        L.b = std::sqrt(L.kbt * L.friction * 0.5 * dt);
      else if (L.kb > 0 && L.temperature > 0)
        L.b = std::sqrt(L.kb * L.temperature * L.friction * 0.5 * dt);
      else
        caviar_->log.error_all(FC_FILE_LINE_FUNC, "expected ('temperature'  & 'kb') or 'kbt' input");
      L.c = 2.0 * dt / (2.0 + L.friction * dt);

      // printing the parameters for langevin

      std::cout << "langevin_parameters:\n"
                << "a: " << L.a << " b: " << L.b << " c: " << L.c << "\n"
                << "b/(0.5*dt): " << L.b / (0.5 * dt) << "(acceleration dimension)\n"
                << "c*a: " << L.c * L.a << "(position dimension)\n"
                << "c*b: " << L.c * L.b << "(position dimension)\n\n"
                << std::flush;

      break;
    }

    initialized = true;
  }

  void Md_simulator::re_calculate_acc()
  {
    // Here there's 'r(t)' stored in 'atom_data->atom_struct_owned.position'
    // First calculate 'a(t)' by using 'r(t)'
    bool update_neighborlist = boundary_condition();

    // for (auto &&n : caviar_->neighborlist)
    //   n->build(update_neighborlist);

    caviar_->neighborlist.build(update_neighborlist);

    atom_data->reset_owned_acceleration();

    for (auto &&f : caviar_->forces)
      f->calculate_acceleration();
  }

  void Md_simulator::step(int64_t i)
  {
    current_step = i;
    step();
  }

  void Md_simulator::step()
  {

    if (!initialized)
    {

      initialize();

      atom_data->reset_virial();

      for (auto &&c : caviar_->constraints)
        c->verify_settings();

      time = dt * initial_step;

      for (auto &&tf : time_function)
        tf->update_time_variable(time);

      for (auto &&tf : time_function_3d)
        tf->update_time_variable(time);

      setup();
    }

    atom_data->reset_virial();

    // atom_data->record_owned_old_data();

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

    /*
    static int my_mpi_rank = comm->me;

    if (my_mpi_rank==0) {

      integrator -> step_part_I ();

      for (auto&& c : constraint)
        c -> step_part_I (i);

      integrator -> step_part_II ();

      for (auto&& c : constraint)
        c -> step_part_II (i);

      bool update_neighborlist = boundary_condition ();

      for (auto &&n : neighborlist)
        n->build(update_neighborlist);

    }

    atom_data -> synch_owned_data(0);

    atom_data -> reset_owned_acceleration();


    for (auto&& f : force_field)
      f -> calculate_acceleration ();


    MPI_Barrier(MPI_COMM_WORLD);

    if (my_mpi_rank==0) {

      integrator -> step_part_III ();

      for (auto&& c : constraint)
        c -> step_part_III (i);
    }

    for (auto&& w : writer)
      w -> write (i, time);
    */

#else

    switch (integrator_type)
    {

    default:
    case Integrator_t::Velocity_verlet:
      integrate_velocity_verlet();
      break;

    case Integrator_t::Leap_frog:
      integrate_leap_frog();
      break;

    case Integrator_t::Velocity_verlet_langevin:
      integrate_velocity_verlet_langevin();
      break;
    }

#endif

    time += dt;
    for (auto &&tf : time_function)
      tf->update_time_variable(time);
    for (auto &&tf : time_function_3d)
      tf->update_time_variable(time);
  }

  bool Md_simulator::boundary_condition()
  {
    bool result = false;

    // the result only in MPI case can
    // become true, due to exchanging a particle between domains
    result = atom_data->exchange_owned(current_step);

    atom_data->exchange_ghost(current_step);
    // #if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)
    //
    //  #elif defined(CAVIAR_WITH_MPI)
    //    MPI_Allreduce(MPI::IN_PLACE, &result, 1, MPI::BOOL, MPI::LOR, MPI_COMM_WORLD);
    //  #endif
    return result;
  }

  void Md_simulator::setup()
  {

    t_start = clock();
    if (use_time)
    {
      initial_step = (int)initial_time / dt;
      final_step = (int)final_time / dt;
    }

    // atom_data->record_owned_old_data();

    // if (neighborlist.size() == 0)
    //   caviar_->log.warning("Md_simulator::setup: neighborlist.size() = 0");
    if (caviar_->forces.size() == 0)
      caviar_->log.warning("Md_simulator::setup: force_field.size() = 0");

#if defined(CAVIAR_SINGLE_MPI_MD_DOMAIN)

    atom_data->synch_owned_data(0);

#endif

    // for (auto &&n : neighborlist)
    //   n->init();
    caviar_->neighborlist.init();

    atom_data->exchange_owned(current_step);

    atom_data->exchange_ghost(current_step);

    // for (auto &&n : neighborlist)
    //   n->build(true);
    caviar_->neighborlist.build(true);

    atom_data->reset_owned_acceleration();

    for (auto &&f : caviar_->forces)
      f->calculate_acceleration();

    atom_data->finalize_temperature();

    // atom_data->finalize_pressure();

    for (auto &&w : caviar_->writers)
      w->write(0, 0.0);
  }

  void Md_simulator::cleanup()
  {
  }

  void Md_simulator::integrate_velocity_verlet()
  {

    auto &pos = atom_data->atom_struct_owned.position;
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &acc = atom_data->atom_struct_owned.acceleration;
    auto &pos_old = atom_data->atom_struct_owned.position_old;
    auto &acc_old = atom_data->atom_struct_owned.acceleration_old; // velocity verlet
    auto psize = pos.size();
    pos_old.resize(psize);
    // auto &pos_old = atom_data -> atom_struct_owned.position_old; //verlet - also shake and m_shake uses this
    // const auto two_dt_inv = 1.0/(2.0*dt);//verlet
    //-------------------------------
    //-------------------------------

    // re_calculate_acc(); // use r(t) to calculate a(t). The forces has to be velocity independent
    //  It has been calculated in previous step

    acc_old.resize(psize);

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < psize; i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      if (std::isnan(acc[i].x) || std::isnan(acc[i].y) || std::isnan(acc[i].z))
      {
        caviar_->log.error_all(FC_FILE_LINE_FUNC,
                               "(MPI_RANK: " + std::to_string(my_mpi_rank) +
                                   ", timestep:" + std::to_string(current_step) + "): The atom with id " + std::to_string(atom_data->atom_struct_owned.id[i]) + " and index " + std::to_string(i) + " at position (" + std::to_string(pos[i].x) +
                                   " , " + std::to_string(pos[i].y) + " , " + std::to_string(pos[i].z) + ") is has NaN acceleration.");
      }
      pos_old[i] = pos[i];
      pos[i] += vel[i] * dt + 0.5 * acc[i] * dt * dt; // r(t+dt) = r(t) + v(t)*dt + 1/2 * a(t) * dt^2
      acc_old[i] = acc[i];                            // calculate before
    }

    for (auto &&c : caviar_->constraints)
    {
      c->fix_position(current_step);
      c->apply_shake(current_step);
    }

    re_calculate_acc(); // use r(t+dt) to calculate a(t+dt)

    for (auto &&c : caviar_->constraints)
      c->fix_acceleration(current_step);

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < vel.size(); i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      vel[i] += 0.5 * dt * (acc[i] + acc_old[i]); // v(t+dt) = v(t) + ( a(t+dt) + a(t) ) * dt / 2
    }

    atom_data->finalize_temperature();

    bool recalculate_temperature = false;

    for (auto &&c : caviar_->constraints)
    {
      c->fix_velocity(current_step, recalculate_temperature);
      c->apply_thermostat(current_step, recalculate_temperature);
    }

    if (recalculate_temperature)
      atom_data->finalize_temperature();

    atom_data->finalize_pressure();

    bool fix_position_needed = false;
    for (auto &&c : caviar_->constraints)
    {
      c->apply_barostat(current_step, fix_position_needed);
    }

    if (fix_position_needed && atom_data->get_pressure_process())
    {
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
      for (unsigned int i = 0; i < psize; i++)
      {
#ifdef CAVIAR_WITH_MPI
        if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
          continue;
#endif
        pos_old[i] = pos[i]; // If not fixed, an artifact velocity will be added to the fixed velocity.
      }

      // double old_virial = atom_data->virialConstraint;

      // for (auto &&c : constraint)
      //   c->apply_shake(current_step); // Shake positions of the atoms again after re-scaling.

      // atom_data->virialConstraint = old_virial; // It is useful if one wants to output it in a file
    }

    for (auto &&c : caviar_->constraints)
      c->apply(current_step);

    for (auto &&w : caviar_->writers)
      w->write(current_step, time); // pos = r(t+dt) , vel = v(t+dt)
  }

  void Md_simulator::integrate_velocity_verlet_langevin()
  {
    auto &pos = atom_data->atom_struct_owned.position;
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &acc = atom_data->atom_struct_owned.acceleration;
    auto &pos_old = atom_data->atom_struct_owned.position_old;
    auto psize = pos.size();
    pos_old.resize(psize);
    // auto &pos_old = atom_data -> atom_struct_owned.position_old; //verlet - also shake and m_shake uses this
    // const auto two_dt_inv = 1.0/(2.0*dt);//verlet

    auto &L = langevin_param;

    L.eta_x.resize(psize);
    L.eta_y.resize(psize);
    L.eta_z.resize(psize);

    for (unsigned int i = 0; i < psize; i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      L.eta_x[i] = L.rnd_ndist_x(L.rnd_tools_x);
      L.eta_y[i] = L.rnd_ndist_y(L.rnd_tools_y);
      L.eta_z[i] = L.rnd_ndist_z(L.rnd_tools_z);
    }

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < psize; i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      if (std::isnan(acc[i].x) || std::isnan(acc[i].y) || std::isnan(acc[i].z))
      {
        caviar_->log.error_all(FC_FILE_LINE_FUNC,
                               "(MPI_RANK: " + std::to_string(my_mpi_rank) +
                                   ", timestep:" + std::to_string(current_step) + "): The atom with id " + std::to_string(atom_data->atom_struct_owned.id[i]) + " and index " + std::to_string(i) + " at position (" + std::to_string(pos[i].x) +
                                   " , " + std::to_string(pos[i].y) + " , " + std::to_string(pos[i].z) + ") is has NaN acceleration.");
      }

      const auto eta = Vector3d<double>{L.eta_x[i], L.eta_y[i], L.eta_z[i]};

      vel[i] += 0.5 * acc[i] * dt + L.b * eta;
      // std::cout << "acc " << acc[i] << ",vel " << vel[i] << "\n";
    }

    atom_data->finalize_temperature();

    bool recalculate_temperature = false;

    for (auto &&c : caviar_->constraints)
    {
      c->fix_velocity(current_step, recalculate_temperature);
      c->apply_thermostat(current_step, recalculate_temperature);
    }

    if (recalculate_temperature)
      atom_data->finalize_temperature();

    atom_data->finalize_pressure();

    bool fix_position_needed = false;
    for (auto &&c : caviar_->constraints)
    {
      c->apply_barostat(current_step, fix_position_needed);
    }

    if (fix_position_needed && atom_data->get_pressure_process())
    {
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
      for (unsigned int i = 0; i < psize; i++)
      {
#ifdef CAVIAR_WITH_MPI
        if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
          continue;
#endif
        pos_old[i] = pos[i]; // If not fixed, an artifact velocity will be added to the fixed velocity.
      }

      // double old_virial = atom_data->virialConstraint;

      // for (auto &&c : constraint)
      //   c->apply_shake(current_step); // Shake positions of the atoms again after re-scaling.

      // atom_data->virialConstraint = old_virial; // It is useful if one wants to output it in a file
    }

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < psize; i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      // std::cout << "b " << pos[i] ;
      pos[i] += vel[i] * L.c;
      // std::cout << "a " << pos[i] << "\n";
    }

    for (auto &&c : caviar_->constraints)
    {
      c->fix_position(current_step);
      c->apply_shake(current_step);
    }

    re_calculate_acc(); // use r(t+dt) to calculate a(t+dt)

    for (auto &&c : caviar_->constraints)
      c->fix_acceleration(current_step);

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < vel.size(); i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      const auto eta = Vector3d<double>{L.eta_x[i], L.eta_y[i], L.eta_z[i]};
      vel[i] = L.a * vel[i] + L.b * eta + 0.5 * acc[i] * dt;
      // std::cout << "2 acc " << acc[i] << ",vel " << vel[i] << "\n";
    }

    for (auto &&c : caviar_->constraints)
      c->apply(current_step);

    for (auto &&w : caviar_->writers)
      w->write(current_step, time); // pos = r(t+dt) , vel = v(t+dt)
  }

  void Md_simulator::integrate_leap_frog()
  {
    auto &pos = atom_data->atom_struct_owned.position;
    auto &pos_old = atom_data->atom_struct_owned.position_old;
    auto &vel = atom_data->atom_struct_owned.velocity;
    auto &acc = atom_data->atom_struct_owned.acceleration;

    // auto &acc_old = atom_data -> atom_struct_owned.acceleration_old; //velocity verlet
    // auto &pos_old = atom_data -> atom_struct_owned.position_old; //verlet - also shake and m_shake uses this
    // const auto two_dt_inv = 1.0/(2.0*dt);//verlet
    //-------------------------------
    //-------------------------------
    // re_calculate_acc(); // use r(t) to calculate a(t). The forces has to be velocity independent
    // It has been calculated in previous step
    auto psize = pos.size();
    pos_old.resize(psize);

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < psize; i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      if (std::isnan(acc[i].x) || std::isnan(acc[i].y) || std::isnan(acc[i].z))
      {
        caviar_->log.error_all(FC_FILE_LINE_FUNC,
                               "(MPI_RANK: " + std::to_string(my_mpi_rank) +
                                   ", timestep:" + std::to_string(current_step) + "): The atom with id " + std::to_string(atom_data->atom_struct_owned.id[i]) + " and index " + std::to_string(i) + " at position (" + std::to_string(pos[i].x) +
                                   " , " + std::to_string(pos[i].y) + " , " + std::to_string(pos[i].z) + ") has NaN acceleration.");
      }
      vel[i] += 0.5 * dt * acc[i]; // v (t + dt/2) = v (t) + (dt/2) a (t)
      pos_old[i] = pos[i];
      pos[i] += dt * vel[i]; // r (t + dt) = r (t) + dt * v (t + dt/2)
    }

    for (auto &&c : caviar_->constraints)
    {
      c->fix_position(current_step);
      c->apply_shake(current_step);
    }

    re_calculate_acc(); // use r(t+dt) to calculate a(t+dt)

    for (auto &&c : caviar_->constraints)
      c->fix_acceleration(current_step);

#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
    for (unsigned int i = 0; i < vel.size(); i++)
    {
#ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
        continue;
#endif
      vel[i] += 0.5 * dt * acc[i]; // v (t + dt) = v (t + dt/2) + (dt/2) a (t + dt)
    }

    atom_data->finalize_temperature();

    bool recalculate_temperature = false;

    for (auto &&c : caviar_->constraints)
    {
      c->fix_velocity(current_step, recalculate_temperature);
      c->apply_thermostat(current_step, recalculate_temperature);
    }

    if (recalculate_temperature)
      atom_data->finalize_temperature();

    atom_data->finalize_pressure();

    bool fix_position_needed = false;
    for (auto &&c : caviar_->constraints)
    {
      c->apply_barostat(current_step, fix_position_needed);
    }

    if (fix_position_needed && atom_data->get_pressure_process())
    {
#ifdef CAVIAR_WITH_OPENMP
#pragma omp parallel for
#endif
      for (unsigned int i = 0; i < psize; i++)
      {
#ifdef CAVIAR_WITH_MPI
        if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank)
          continue;
#endif
        pos_old[i] = pos[i]; // If not fixed, an artifact velocity will be added to the fixed velocity.
      }

      // double old_virial = atom_data->virialConstraint;

      // for (auto &&c : constraint)
      //   c->apply_shake(current_step); // Shake positions of the atoms again after re-scaling.

      // atom_data->virialConstraint = old_virial; // It is useful if one wants to output it in a file
    }

    for (auto &&c : caviar_->constraints)
      c->apply(current_step);

    for (auto &&w : caviar_->writers)
      w->write(current_step, time); // pos = r(t+dt) , vel = v(t+dt)
  }

  /*
  void Md_simulator::integrate_verlet () {
  {
    caviar_->log.error_all(FILE_LINE_FUNC,"not implemented");

    if (pos_old.size() != pos.size()) pos_old.resize(pos.size());

    atom_data->record_owned_position_old = true;

    re_calculate_acc(); // use r(t) to calculate a(t). The forces has to be velocity independent


  #ifdef CAVIAR_WITH_OPENMP
    #pragma omp parallel for
  #endif
    for (unsigned int i=0; i<psize; i++) {
      #ifdef CAVIAR_WITH_MPI
      if (atom_data->atom_struct_owned.mpi_rank[i] != my_mpi_rank) continue;
      #endif
      pos [i] = (2.0*pos[i]) - pos_old[i]  + acc [i] * dt * dt;  // r (t + dt) = 2*r(t) − r (t − dt) + a(t) * dt^2


      // XXX note that this is 'v(t)' not 'v(t+dt)'
      vel [i] = domain->fix_distance(pos[i] - pos_old[i]) * two_dt_inv; // v(t) = (r (t + dt) - r (t - dt) )/ (2*dt)

    }

    // velocity related constraints are useless here
    for (auto&& c : constraint) c -> apply (i); // use unconstrained r(t+dt) and change it to constrained value


    for (auto&& w : writer) w -> write (i, time); // pos = r(t+dt) , vel = v(t)

  }*/
}
