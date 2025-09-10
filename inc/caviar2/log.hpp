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

#include <string>


namespace caviar2
{

#define FC_FILE_LINE_FUNC_LINE_COL __FILE__, __LINE__, __func__, line, col
#define FC_FILE_LINE_FUNC_PARSE __FILE__, __LINE__, __func__, parser->line, parser->col
#define FC_FILE_LINE_FUNC __FILE__, __LINE__, __func__

  /**
   * This class does all of the output massages.
   *
   *
   */
  class Log
  {
  public:
    Log(class Caviar2 *caviar) ;
    Log();
    virtual void set_caviar(Caviar2 *caviar);
    void comment(const std::string &, const bool endline = true);
    void comment(const char *, const bool endline = true);

    void info(const std::string &, const int level = 0, const bool endline = true);
    void info(const char *, const int level = 0, const bool endline = true);

    void info_create(const std::string &, const int level = 1, const bool endline = true);
    void info_create(const char *, const int level = 1, const bool endline = true);

    void info_read(const std::string &, const int level = 1, const bool endline = true);
    void info_read(const char *, const int level = 1, const bool endline = true);

    void warning(const std::string &, const int level = 0, const bool endline = true);
    void warning(const char *, const int level = 0, const bool endline = true);

    void error_all(const std::string &);

    void error_all(const char *, int, const char *, const std::string &, unsigned int, const char *);
    void error_one(const char *, int, const char *, const std::string &, unsigned int, const char *);

    void error_all(const char *, int, const char *, const std::string &, unsigned int, const std::string &);
    void error_one(const char *, int, const char *, const std::string &, unsigned int, const std::string &);

    void error_all(const char *, int, const char *, const std::string &);
    void error_one(const char *, int, const char *, const std::string &);

    bool output_info[5], output_warning[5];

    bool err_flag = true;
    bool log_flag = true;

    Caviar2 *caviar_ = nullptr;
  };

}
