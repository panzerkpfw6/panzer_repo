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
/// @file src/source.c
/// @brief This file contains the implementation of the source simulator.
///
#include <math.h>
#include <stdio.h>
#include <stencil/config.h>
#include <stencil/macros.h>
#include <stencil/shot.h>

void shot_init(shot_t *shot, bool cpu, bool modeling) {
    MSG("OUTDIR=%s",OUTDIR);
    char tmp[512];
    MSG("breakpoint 0");
    sprintf(tmp, "mkdir -p %s", OUTDIR);
    system(tmp);
    sprintf(tmp, "%s/%s_%d.raw", OUTDIR, SNAP_BASE, shot->id);
    MSG("tmp in shot_init=%s",tmp);
    if (!modeling) {
        shot->fd_snap = fopen(tmp, "wb+");
        CHK(shot->fd_snap == NULL, "failed to open snapshot file");
    } else {
        shot->fd_snap = NULL;
    }

#ifdef __DEBUG
    shot->fd_fwd = NULL;
  shot->fd_bwd = NULL;

  if (cpu)
    sprintf(tmp, "%s/%s_fwd_%d.raw", OUTDIR, SNAP_BASE, shot->id);
  else
    sprintf(tmp, "%s/gpu_%s_fwd_%d.raw", OUTDIR, SNAP_BASE, shot->id);
  shot->fd_fwd = fopen(tmp, "wb");
  CHK(shot->fd_fwd == NULL, "fwd snapshot file impossible to open");
  if (!modeling) {
    if (cpu)
      sprintf(tmp, "%s/%s_bwd_%d.raw", OUTDIR, SNAP_BASE, shot->id);
    else
      sprintf(tmp, "%s/gpu_%s_bwd_%d.raw", OUTDIR, SNAP_BASE, shot->id);
    shot->fd_bwd = fopen(tmp, "wb");
    CHK(shot->fd_bwd == NULL, "bwd snapshot file impossible to open");
  }
#endif // __DEBUG

    shot->fd_img = NULL;
    shot->fd_ilm = NULL;
    if (!modeling) {
        if (cpu)
            sprintf(tmp, "%s/%s_%d.raw", OUTDIR, IMG_BASE, shot->id);
        else
            sprintf(tmp, "%s/gpu_%s_%d.raw", OUTDIR, IMG_BASE, shot->id);
        shot->fd_img = fopen(tmp, "wb");
        CHK(shot->fd_img == NULL, "img file impossible to open, aborting");
        if (cpu)
            sprintf(tmp, "%s/%s_%d.raw", OUTDIR, ILM_BASE, shot->id);
        else
            sprintf(tmp, "%s/gpu_%s_%d.raw", OUTDIR, ILM_BASE, shot->id);
        shot->fd_ilm = fopen(tmp, "wb");
        CHK(shot->fd_ilm == NULL, "ilm file impossible to open, aborting");
    }
    MSG("... dealing with shot number %u (source @ %u).", shot->id,shot->srcidx);
}

void shot_release(shot_t *shot) {
  if (shot->fd_snap) {
    char tmp[512];
    sprintf(tmp, "%s/%s_%d.raw", OUTDIR, SNAP_BASE, shot->id);
    MSG("shot release=%s\n",tmp);
    fclose(shot->fd_snap); shot->fd_snap = NULL;
    int err = remove(tmp);
    CHK(err != 0, "failed to delete snapshot file, aborting");
  }
//  char tmp[512];
//  sprintf(tmp, "%s/%s_%d.raw", OUTDIR, SNAP_BASE, shot->id);
//  MSG("shot release=%s\n",tmp);

//  MSG("shot->fd_img=%s\n",shot->fd_img);
//  MSG("shot->fd_ilm=%s\n",shot->fd_ilm);

  if (shot->fd_img) { fclose(shot->fd_img); shot->fd_img = NULL; }
  if (shot->fd_ilm) { fclose(shot->fd_ilm); shot->fd_ilm = NULL; }

  #ifdef __DEBUG
  if (shot->fd_fwd) fclose(shot->fd_fwd);
  if (shot->fd_bwd) fclose(shot->fd_bwd);
  #endif // __DEBUG
}
