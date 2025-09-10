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

#include "caviar2/log.hpp"
#include <iostream>
#include <fstream>

#if defined(CAVIAR_WITH_MPI)
#include <mpi.h>
#endif
namespace caviar2
{

  Log::Log(Caviar2 *fptr)
  {
    caviar_ = fptr;
    set_caviar(fptr);
    for (int i = 0; i < 5; ++i)
    {
      output_info[i] = true;
      output_warning[i] = true;
    }
  }

  Log::Log()
  {

    for (int i = 0; i < 5; ++i)
    {
      output_info[i] = true;
      output_warning[i] = true;
    }
  };

  void Log::set_caviar(Caviar2 *c)
  {
  }

  void Log::info(const char *str, int level, bool endline)
  {
    info(static_cast<std::string>(str), level, endline);
  }
  void Log::info(const std::string &str, int level, bool endline)
  {
    if (output_info[level])
    {
#if defined(CAVIAR_WITH_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
      int me = comm->me;
      if (me == 0)
      {
#endif
        std::cout << "[INF] ";
        std::cout << str;
        if (endline)
          std::cout << std::endl;
        else
          std::cout << std::flush;
#if defined(CAVIAR_WITH_MPI)
      }
#endif
    }
  }

  void Log::info_create(const char *str, int level, bool endline)
  {
    info_create(static_cast<std::string>(str), level, endline);
  }
  void Log::info_create(const std::string &str, int level, bool endline)
  {
    std::string s = "(Create) " + str;
    info(s, level, endline);
  }

  void Log::info_read(const char *str, int level, bool endline)
  {
    info_read(static_cast<std::string>(str), level, endline);
  }
  void Log::info_read(const std::string &str, int level, bool endline)
  {
    std::string s = "(Call) " + str + ".read()";
    info(s, level, endline);
  }

  void Log::warning(const char *str, int level, bool endline)
  {
    warning(static_cast<std::string>(str), level, endline);
  }
  void Log::warning(const std::string &str, int level, bool endline)
  {
    if (output_warning[level])
    {
#if defined(CAVIAR_WITH_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
      int me = comm->me;
      if (me == 0)
      {
#endif
        std::cout << "[WRN] ";
        std::cout << str;
        if (endline)
          std::cout << std::endl;
        else
          std::cout << std::flush;
#if defined(CAVIAR_WITH_MPI)
      }
#endif
    }
  }

  void Log::comment(const char *str, bool endline)
  {
    comment(static_cast<std::string>(str), endline);
  }
  void Log::comment(const std::string &str, bool endline)
  {
#if defined(CAVIAR_WITH_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    int me = comm->me;
    if (me == 0)
    {
#endif
      std::cout << str;
      if (endline)
        std::cout << std::endl;
      else
        std::cout << std::flush;
#if defined(CAVIAR_WITH_MPI)
    }
#endif
  }

  /*
    bool Log::read(caviar2::interpreter::Parser *parser)
    {
      info("output read");
      bool in_file = true;
      while (true)
      {
        GET_A_TOKEN_FOR_CREATION
        auto t = token.string_value;
        if (string_cmp(t, "info"))
        {
          int level = -1;
          GET_OR_CHOOSE_A_INT(level, "", "")
          bool status = 0;
          GET_A_BOOL(status, "", "")
          if (level > 5 || level < 0)
            caviar_->log.error_all(FC_FILE_LINE_FUNC_PARSE, "'info' level is defined between 0 to 4. To choose all, use 5.");
          if (level == 5)
            for (int i = 0; i < 5; ++i)
              output_info[i] = status;
          else
            output_info[level] = status;
        }
        else if (string_cmp(t, "warning"))
        {
          int level = -1;
          GET_OR_CHOOSE_A_INT(level, "", "")
          bool status = 0;
          GET_A_BOOL(status, "", "")
          if (level > 5 || level < 0)
            caviar_->log.error_all(FC_FILE_LINE_FUNC_PARSE, "'warning' level is defined between 0 to 4. To choose all, use 5.");
          if (level == 5)
            for (int i = 0; i < 5; ++i)
              output_warning[i] = status;
          else
            output_warning[level] = status;
        }
        else
          caviar_->log.error_all(FC_FILE_LINE_FUNC_PARSE, "Invalid syntax: This output command doesn't exist.");
      }
      return in_file;


    }*/

  // All procs must call this else there would be a deadlock
  void Log::error_all(const std::string &str)
  {
    std::cerr << "[ERR] " << str << std::endl;
    exit(1);
  }

  void Log::error_all(const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const char *str)
  {
    int me = 0;
#ifdef CAVIAR_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (me == 0)
    {
      if (err_flag)
      {
        std::cerr << "[ERR] " << str << std::endl;
        std::cerr << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
        std::cerr << parsing_line << std::endl;
        for (unsigned int i = 0; i < col; ++i)
          std::cerr << ' ';
        std::cerr << '^' << std::endl;
      }
      if (log_flag)
      {
        std::cout << "[ERR] " << str << std::endl;
        std::cout << " '" << file << ':' << line << " in '" << func << "'." << std::endl;
        std::cout << parsing_line << std::endl;
        for (unsigned int i = 0; i < col; ++i)
          std::cout << ' ';
        std::cout << '^' << std::endl;
      }
      // if (log_flag)
      //   log.close();
    }
#ifdef CAVIAR_WITH_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  // One proc calling this will abort all

  void Log::error_one(const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const char *str)
  {
#ifdef CAVIAR_WITH_MPI
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (err_flag)
    {
      std::cerr << "[ERR] " << str << std::endl;
      std::cerr << " MPI rank " << me << std::endl;
      std::cerr << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;

      std::cerr << parsing_line << std::endl;
      for (unsigned i = 0; i < col; ++i)
        std::cerr << ' ';
      std::cerr << '^' << std::endl;
    }
    if (log_flag)
    {
      std::cout << "[ERR] " << str << std::endl;
      std::cerr << " MPI rank " << me << std::endl;
      std::cout << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      std::cout << parsing_line << std::endl;
      for (unsigned i = 0; i < col; ++i)
        std::cout << ' ';
      std::cout << '^' << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    error_all(file, line, func, parsing_line, col, str);
#endif
  }

  void Log::error_all(const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const std::string &str)
  {
    int me = 0;
#ifdef CAVIAR_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    //  /*
    if (me == 0)
    {
      if (err_flag)
      {
        std::cerr << "[ERR] " << str << std::endl;
        std::cerr << " MPI rank " << me << std::endl;
        std::cerr << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
        std::cerr << parsing_line << std::endl;
        for (unsigned int i = 0; i < col; ++i)
          std::cerr << ' ';
        std::cerr << '^' << std::endl;
      }
      if (log_flag)
      {
        std::cout << "[ERR] " << str << std::endl;
        std::cerr << " MPI rank " << me << std::endl;
        std::cout << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
        std::cout << parsing_line << std::endl;
        for (unsigned int i = 0; i < col; ++i)
          std::cout << ' ';
        std::cout << '^' << std::endl;
      }
      // if (log_flag)
      //   log.close();
    }
//  */
#ifdef CAVIAR_WITH_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  void Log::error_one(const char *file, int line, const char *func, const std::string &parsing_line, unsigned int col, const std::string &str)
  {
#ifdef CAVIAR_WITH_MPI
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (err_flag)
    {
      std::cerr << "[ERR] " << str << std::endl;
      std::cerr << " MPI rank " << me << std::endl;
      std::cerr << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      std::cerr << parsing_line << std::endl;
      for (unsigned i = 0; i < col; ++i)
        std::cerr << ' ';
      std::cerr << '^' << std::endl;
    }
    if (log_flag)
    {
      std::cout << "[ERR] " << str << std::endl;
      std::cerr << " MPI rank " << me << std::endl;
      std::cout << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      std::cout << parsing_line << std::endl;
      for (unsigned i = 0; i < col; ++i)
        std::cout << ' ';
      std::cout << '^' << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    error_all(file, line, func, parsing_line, col, str);
#endif
  }

  void Log::error_all(const char *file, int line, const char *func, const std::string &str)
  {
    int me = 0;
#ifdef CAVIAR_WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (me == 0)
    {
      if (err_flag)
      {
        std::cerr << "[ERR] " << str << std::endl;
        std::cerr << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      }
      if (log_flag)
      {
        std::cout << "[ERR] " << str << std::endl;
        std::cout << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
      }
      // if (log_flag)
      //   log.close();
    }

#ifdef CAVIAR_WITH_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  void Log::error_one(const char *file, int line, const char *func, const std::string &str)
  {
#ifdef CAVIAR_WITH_MPI
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (err_flag)
    {
      std::cerr << "[ERR] " << str << std::endl;
      std::cerr << " MPI rank " << me << std::endl;
      std::cerr << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
    }
    if (log_flag)
    {
      std::cout << "[ERR] " << str << std::endl;
      std::cout << " MPI rank " << me << std::endl;
      std::cout << " '" << file << ':' << line << "' in '" << func << "'." << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
#else
    error_all(file, line, func, str);
#endif
  }
}
