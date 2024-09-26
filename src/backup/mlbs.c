#include <semaphore.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include <stencil/macros.h>
#include <stencil/wave_tb.h>

typedef struct mlbs_s mlbs_t;

struct mlbs_s {
  sem_t * task;
  pthread_mutex_t mutex;

  int nb_total_task;
  int count;

  int y_len_r;
  int y_len_l;

  float ** ldata;
  float ** rdata;

  sizt_t datasize;

  volatile int * index; // 1000 is left, 2000 is right
  volatile int head;
  volatile int tail;
  int indexsize;

  // writer
  size_t * ldataoffset;
  size_t * rdataoffset;


  // reader
};

int mlbs_init(mlbs_t * ctx,
              tb_t * tb) {
  int rc;

  rc = sem_init(ctx->task,0,0);
  if (rc != 0) ERR_MSG("sem_init %s\n",strerror(errno));

  ctx->ldata = calloc(ctx->y_len_l,sizeof(float*));
  ctx->rdata = calloc(ctx->y_len_r,sizeof(float*));

  ctx->ldataoffset = calloc(ctx->y_len_l,sizeof(size_t));
  ctx->rdataoffset = calloc(ctx->y_len_r,sizeof(size_t));

  size_t fwdsize = 1LL * ctx->nnx * ctx->nny * ctx->diam_width;
  ctx->datasize = fwdsize;

  for (int i = 0; i < ctx->y_len_l; i ++) {
    ctx->ldata[i] = calloc(fwdsize,sizeof(float));
  }

  for (int i = 0; i < ctx->y_len_r - 2; i ++) {
    ctx->rdata[i] = calloc(fwdsize,sizeof(float));
  }

  ctx->rdata[ctx->y_len_r - 2] = calloc(fwdsize/2,sizeof(float));
  ctx->rdata[ctx->y_len_r - 1] = calloc(fwdsize/2,sizeof(float));

  ctx->nb_total_task = ; // TODO

  ctx->indexsize = ctx->y_len_l + ctx->y_len_r;
  ctx->index = calloc(ctx->indexsize,sizeof(int));
  ctx->head = 0;
  ctx->tail = 0;

  rc = pthread_mutex_init(&ctx->mutex,NULL);
  if (rc != 0) ERR_MSG("pthread_mutex_init %s\n",strerror(errno));

  ctx->y_len_r = tb->y_len_r;
  ctx->y_len_l = tb->y_len_l;
}

void * mlbs_writer_helper(void * data) {
  int rc,index;
  mlbs_t * ctx = (mlbs_t *) data;
  size_t size_done;

  while (data->count != data->nb_total_task) {

    rc = sem_wait(ctx->task);
    if (rc != 0) ERR_MSG("sem_wait %s\n",strerror(errno));

    rc = pthread_mutex_lock(&ctx->mutex);
    if (rc != 0) ERR_MSG("pthread_mutex_lock %s\n",strerror(errno));

    index = ctx->index[ctx->head];
    ctx->head ++;

    rc = pthread_mutex_unlock(&ctx->mutex);
    if (rc != 0) ERR_MSG("pthread_mutex_unlock %s\n",strerror(errno));

    if (index >= 1000 && index < 2000) { // left
      index -= 1000;

      rc = fseek(fp,ctx->ldataoffset[index], SEEK_SET);
      if (rc != 0) ERR_MSG("fseek %s\n",strerror(errno));

      size_done = fwrite(ctx->ldata[index],sizeof(float),ctx->datasize,ctx->fp);
      if (size_done != ctx->datasize) ERR_MSG("fwrite %s\n",strerror(errno));
    } else { // right
      index -= 2000;

      rc = fseek(fp,ctx->rdataoffset[index], SEEK_SET);
      if (rc != 0) ERR_MSG("fseek %s\n",strerror(errno));

      size_t datasize = ctx->datasize;
      if (index >= ctx->y_len_r-2) datasize /= 2;

      rc = fseek(fp,ctx->rdataoffset[index], SEEK_SET);
      if (rc != 0) ERR_MSG("fseek %s\n",strerror(errno));

      size_done = fwrite(ctx->rdata[index],sizeof(float),datasize,ctx->fp);
      if (size_done != ctx->datasize) ERR_MSG("fwrite %s\n",strerror(errno));
    }

    data->count ++;
  }


}

int mlbs_free(mlbs_t * ctx,
              tb_t * tb) {


  rc = sem_destroy(ctx->task);
  if (rc != 0) ERR_MSG("sem_destroy %s\n",strerror(errno));

  for (int i = 0; i < ctx->y_len_l; i ++) {
    free(ctx->ldata[i]);
  }

  for (int i = 0; i < ctx->y_len_r; i ++) {
    free(ctx->rdata[i]);
  }

  free(ctx->rdata); ctx->rdata = NULL;
  free(ctx->ldata); ctx->ldata = NULL;

  free(ctx->index); ctx->index = NULL;

  rc = pthread_mutex_destroy(&ctx->mutex,NULL);
  if (rc != 0) ERR_MSG("pthread_mutex_destroy %s\n",strerror(errno));

  free(ctx->rdataoffset); ctx->rdataoffset = NULL;
  free(ctx->ldataoffset); ctx->ldataoffset = NULL;

}