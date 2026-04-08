#ifndef __STENCIL_MACROS_H_
#define __STENCIL_MACROS_H_
///
/// @copyright Copyright 2024- Pavel Plotnitskii. All rights reserved.
/// This file is part of the \b stencil project.
///
/// \b stencil is free software: you may redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// The stencil project is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with the \b stencil project. If not, see <http://www.gnu.org/licenses/>.
///
/// @author Pavel Plotnitskii
/// @file stencil/macros.h
/// @brief Contains the major settings of project.
///
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>


/// @brief print a formatted message.
#define MSG(fmt,...)                    \
  do {                                  \
    fprintf(stdout, "[STENCIL MSG]:");  \
    fprintf(stdout, fmt, ##__VA_ARGS__);\
    fprintf(stdout, "\n");              \
  } while(0)

/// @brief print a formatted message without line return.
#define MSG_NLR(fmt,...) \
  do {                                  \
    fprintf(stdout, "[STENCIL MSG]:");  \
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
    	      "[STENCIL ERR] %s @ %s:%d\n", \
            str, __FILE__, __LINE__);     \
    exit(EXIT_FAILURE);                   \
  }                                       \
} while(0)

/// @brief print an error message and fail.
#define ERR(str)                        \
do {                                    \
  fprintf(stderr,                       \
    	    "[STENCIL ERR] %s @ %s:%d\n", \
          str, __FILE__, __LINE__);     \
  exit(EXIT_FAILURE);                   \
} while(0)

/// @brief print an error formatted message.
#define ERR_MSG(fmt,...)                     \
  do {                                       \
    fprintf(stderr, "[STENCIL ERR @ %s:%d]", \
           __FILE__,__LINE__ );              \
    fprintf(stderr, fmt, ##__VA_ARGS__);     \
    exit(EXIT_FAILURE);                      \
  } while(0)

/// @brief print an error message and fail if predicate.
#define ERR_IF(predicate, str)            \
do {                                      \
  if(predicate) {                         \
    fprintf(stderr,                       \
    	      "[STENCIL ERR] %s @ %s:%d\n", \
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
            "[STENCIL ERR] failed to create buffer @ %s:%d\n",    \
            __FILE__, __LINE__);                                  \
    exit(EXIT_FAILURE);                                           \
  }                                                               \
  memset(buffer, 0, (size)*sizeof(float))

  #define CREATE_BUFFER_ONLY(buffer, size)                               \
  CHK(posix_memalign((void**)&buffer, 4096, (size)*sizeof(float)),\
      "failed to allocate heap memory");                          \
  if (buffer == NULL) {                                           \
    fprintf(stderr,                                               \
            "[STENCIL ERR] failed to create buffer @ %s:%d\n",    \
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
    fprintf(stderr, "[STENCIL GPU ERR] %s @ %s:%d\n",             \
            cudaGetErrorString(status),      __FILE__, __LINE__); \
    exit(EXIT_FAILURE);                                           \
  }                                                               \
} while(0);

//////////  special TB macros   //////////
// define type
typedef int      Myint ;
typedef float    Myfloat ;
typedef long int Myint64 ;
typedef float hFloat ;
typedef float real_t;

//#define U1(i,j,k)         (u1[((1ULL)*((i)*(nny)+(j))*(nnz)+(k))])
//#define U2(i,j,k)         (u2[((1ULL)*((i)*(nny)+(j))*(nnz)+(k))])
//#define U3(i,j,k)         (u3[((1ULL)*((i)*(nny)+(j))*(nnz)+(k))])
#define V1_xyz(i,j,k)         (v1[((1ULL)*((i)*(nny)+(j))*(nnz)+(k))])
//#define V2(i,j,k)         (v2[((1ULL)*((i)*(nny)+(j))*(nnz)+(k))])
//#define V3(i,j,k)         (v3[((1ULL)*((i)*(nny)+(j))*(nnz)+(k))])
//#define ROC(i,j,k)        (roc2[((1ULL)*((i)*(nny)+(j))*(nnz)+(k))])
//////////////////////////////
#define MAX(a,b)    (a>b?a:b)
#define MIN(a,b)    (a<b?a:b)
#define max(a, b) ((a) > (b) ? (a) : (b))

#endif //  __STENCIL_MACROS_H_
