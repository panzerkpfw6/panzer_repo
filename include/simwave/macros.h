#ifndef __SIMWAVE_MACROS_H_
#define __SIMWAVE_MACROS_H_
///
/// @copyright Copyright 2017- Issam Said. All rights reserved.
/// This file is part of \b simwave.
///
/// \b simwave is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// simwave is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with \b simwave.  If not, see <http://www.gnu.org/licenses/>.
///
/// @author Issam Said
/// @file simwave/macros.h
/// @brief Contains the major settings of project.
///
#include <stdlib.h>

/// @brief print a formatted message.
#define MSG(fmt,...)                    \
  do {                                  \
    fprintf(stdout, "[SIMWAVE MSG]:");  \
    fprintf(stdout, fmt, ##__VA_ARGS__);\
    fprintf(stdout, "\n");              \
  } while(0)

/// @brief print a formatted message without line return.
#define MSG_NLR(fmt,...) \
  do {                                  \
    fprintf(stdout, "[SIMWAVE MSG]:");  \
    fprintf(stdout, fmt, ##__VA_ARGS__);\
  } while(0)

/// @brief print a non-formatted message without line return.
#define MSG_INL(fmt,...) fprintf(stdout, fmt, ##__VA_ARGS__)

/// @brief print a non-formatted line return.
#define MSG_JLR() fprintf(stdout, "\n")

/// @brief check the return code of a given call and fail if error.
#define CHK(call, str)                    \
do {                                      \
  if(call) {                              \
    fprintf(stderr,                       \
    	      "[SIMWAVE ERR] %s @ %s:%d\n", \
            str, __FILE__, __LINE__);     \
    exit(EXIT_FAILURE);                   \
  }                                       \
} while(0)

/// @brief print an error message and fail.
#define ERR(str)                        \
do {                                    \
  fprintf(stderr,                       \
    	    "[SIMWAVE ERR] %s @ %s:%d\n", \
          str, __FILE__, __LINE__);     \
  exit(EXIT_FAILURE);                   \
} while(0)

/// @brief print an error formatted message.
#define ERR_MSG(fmt,...)                     \
  do {                                       \
    fprintf(stderr, "[SIMWAVE ERR @ %s:%d]", \
           __FILE__,__LINE__ );              \
    fprintf(stderr, fmt, ##__VA_ARGS__);     \
    exit(EXIT_FAILURE);                      \
  } while(0)

/// @brief print an error message and fail if predicate.
#define ERR_IF(predicate, str)            \
do {                                      \
  if(predicate) {                         \
    fprintf(stderr,                       \
    	      "[SIMWAVE ERR] %s @ %s:%d\n", \
            str, __FILE__, __LINE__);     \
    exit(EXIT_FAILURE);                   \
  }                                       \
} while(0)

#define MAX(a,b)    (a>b?a:b)
#define MIN(a,b)    (a<b?a:b)
#define CREATE_BUFFER(buffer, size)                               \
  CHK(posix_memalign((void**)&buffer, 4096, (size)*sizeof(float)),\
      "failed to allocate heap memory");                          \
  if (buffer == NULL) {                                           \
    fprintf(stderr,                                               \
            "[SIMWAVE ERR] failed to create buffer @ %s:%d\n",    \
            __FILE__, __LINE__);                                  \
    exit(EXIT_FAILURE);                                           \
  }                                                               \
  memset(buffer, 0, (size)*sizeof(float))

  #define CREATE_BUFFER_ONLY(buffer, size)                               \
  CHK(posix_memalign((void**)&buffer, 4096, (size)*sizeof(float)),\
      "failed to allocate heap memory");                          \
  if (buffer == NULL) {                                           \
    fprintf(stderr,                                               \
            "[SIMWAVE ERR] failed to create buffer @ %s:%d\n",    \
            __FILE__, __LINE__);                                  \
    exit(EXIT_FAILURE);                                           \
  }

#define NULIFY_BUFFER(buffer, size) memset(buffer, 0, size*sizeof(float))

#define DELETE_BUFFER(buffer) free(buffer)

/// @brief Checks the returned codes of GPU calls.
#define GPU_CHK(call)                                             \
do {                                                              \
  cudaError_t status = call;                                      \
  if (status != cudaSuccess) {                                    \
    fprintf(stderr, "[SIMWAVE GPU ERR] %s @ %s:%d\n",             \
            cudaGetErrorString(status),      __FILE__, __LINE__); \
    exit(EXIT_FAILURE);                                           \
  }                                                               \
} while(0);

#endif //  __SIMWAVE_MACROS_H_
