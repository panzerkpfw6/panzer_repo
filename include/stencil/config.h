#ifndef __STENCIL_CONFIG_H_
#define __STENCIL_CONFIG_H_
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
/// @file stencil/config.h
/// @brief contains the major settings of project.
///

#ifdef _OPENMP
#include <omp.h>
#endif

/// @def __GETTIMEOFDAY
/// @brief set Linux gettimeofday wall-clock timer
/// to be used by the timing routines \see timer.h.
#define __GETTIMEOFDAY 1

#if defined(__APPLE__) || defined(__MACOSX)
#define __GETTIMEOFDAY 1
#endif

#ifdef __DEBUG

/// @brief dump the source function in a file.
#define __DUMP_SOURCE

/// @brief dump the generated velocity in a file.
#define __DUMP_VEL

/// @brief dump the generated pml in a file.
#define __DUMP_PML

#endif // __DEBUG

/// @def __SOURCE_FILE
/// @brief path to the source file.
#ifdef __DUMP_SOURCE
#define SOURCE_BASE "source"
#endif // __DUMP_SOURCE

#define OUTDIR       "data"
//#define OUTDIR       "/raid/stencil"
#define SNAP_BASE    "snap"
#define SISMOS_BASE  "sismos"
#define VEL          "vel.raw"
#define IMG_BASE     "img"
#define ILM_BASE     "ilm"

#ifdef __DEBUG
#define SNAP_DBG     "snap_dbg"
#endif // __DEBUG

#endif //  __STENCIL_CONFIG_H_
