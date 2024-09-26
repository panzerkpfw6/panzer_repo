#ifndef __STENCIL_SHOT_H_
#define __STENCIL_SHOT_H_

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
/// @file stencil/shot.h
/// @brief A shot descriptor with its manipulation routines.
///
/// Define a structure for the shots and their manipulation routines.
///
#include <stdio.h>
#include <stdbool.h>
#include <stencil/config.h>

typedef struct __shot_t {
	unsigned int  srcidx;
	unsigned int      id;
	FILE *fd_snap;
//	#ifdef __DEBUG
	FILE* fd_fwd;
	FILE* fd_bwd;
	FILE* fd_img;
	FILE* fd_ilm;
//	#endif // __DEBUG
} shot_t;

void shot_init(shot_t *, bool, bool);

void shot_release(shot_t *);

#endif // __STENCIL_SHOT_H_
