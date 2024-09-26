#ifndef __STENCIL_BWRITER_H_
#define __STENCIL_BWRITER_H_

typedef struct bwriter_s bwriter_t;

int bwriter_create(bwriter_t ** ctx);

int bwriter_destroy(bwriter_t ** ctx);

int bwriter_init(bwriter_t * ctx,
                 int nnx,
                 int nny,
                 int diam_width,
                 int num_thread_groups,
                 char * filename);

int bwriter_write(bwriter_t * ctx,
                  int groupid,
                  float  * data,
                  size_t dataoffset,
                  size_t datasize);

int bwriter_free(bwriter_t * ctx);
#endif
