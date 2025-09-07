/*
 * -----------------------------------------------------------------------------
 * CAVIAR2 - C++ Library
 * 
 * Copyright (c) 2025 Morad Biagooi and Ehsan Nedaaee Oskoee
 * All rights reserved.
 * 
 * This software is provided "as is", without warranty of any kind.
 * You may use, copy, modify, and distribute this software for any purpose,
 * provided this copyright notice is retained.
 * 
 * -----------------------------------------------------------------------------
 */

#pragma once

#include <string>

namespace caviar2 {

/**
 * @brief A sample class for caviar2 library.
 * 
 * This class demonstrates a basic interface.
 */
class Caviar2 {
public:
    Caviar2();
    ~Caviar2();

    /**
     * @brief Initialize the library or object.
     * @return true on success, false otherwise.
     */
    bool initialize();

    /**
     * @brief Perform some computation.
     * @param input An input string.
     * @return Result string after processing.
     */
    std::string process(const std::string& input) const;

private:
    // Private implementation details (PIMPL or just members)
    int internal_state_;
};

} // namespace caviar2
