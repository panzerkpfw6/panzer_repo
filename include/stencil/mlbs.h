#ifndef __STENCIL_MLBS_H_
#define __STENCIL_MLBS_H_

typedef struct mlbs_s mlbs_t;

int mlbs_create(mlbs_t ** ctx);

int mlbs_destroy(mlbs_t ** ctx);

int mlbs_init(mlbs_t * ctx,
              int nnx,
              int nny,
              int diam_width,
              int num_thread_groups,
              char * filename,
              int coreid);

int mlbs_write(mlbs_t * ctx,
               int groupid,
               float  * data,
               size_t dataoffset,
               size_t datasize);

int mlbs_free(mlbs_t * ctx);
#endif