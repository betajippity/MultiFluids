// Written by Peter Kutz.
// Shorthand cout for multi-threaded programs.
// The fflush is necessary in Qt Creator.

#pragma once

#include <iostream>
#include <sstream>

#define PRINT(X)					\
do {								\
std::ostringstream stream;			\
stream << X;						\
std::cout << stream.str();			\
fflush(stdout);                     \
} while (0)

#define PRINT_LINE(X)				\
do {								\
std::ostringstream stream;			\
stream << X << std::endl;			\
std::cout << stream.str();			\
fflush(stdout);                     \
} while (0)
