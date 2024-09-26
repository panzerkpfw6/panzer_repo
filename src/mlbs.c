#define _GNU_SOURCE             /* See feature_test_macros(7) */
#include <semaphore.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <pthread.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include <stencil/macros.h>
#include <stencil/wave_tb.h>
#include <stencil/mlbs.h>

struct mlbs_s {
  FILE * fp;
  int coreid;
  pthread_t thread_id;
  sem_t  task;
  pthread_mutex_t mutex;

  pthread_mutex_t * mutex_p;
  pthread_mutex_t * mutex_c;
  pthread_cond_t  * mutex_cond;

  int nb_group;

  float ** data;


  int * inuse;
  int * index; // groupid
  volatile int head;
  volatile int tail;

  int num_thread_groups;

  // writer
  size_t * dataoffset;
  size_t * datasize;

  // reader
};

void * mlbs_writer_helper(void * data);

int mlbs_create(mlbs_t ** ctx) {
  *ctx = calloc(1,sizeof(mlbs_t));
  if (*ctx == NULL) ERR_MSG("calloc\n");

  return 0;
}

int mlbs_destroy(mlbs_t ** ctx) {
  free(*ctx);
  *ctx = NULL;

  return 0;
}

int mlbs_init(mlbs_t * ctx,
              int nnx,
              int nny,
              int diam_width,
              int num_thread_groups,
              char * filename,
              int coreid) {
  int rc;

  ctx->fp = fopen(filename,"wb+");
  if (ctx->fp == NULL) ERR_MSG("fopen %s\n",strerror(errno));

  // Remove system buffer for fwrite
  //rc = setvbuf(ctx->fp, NULL, _IONBF,0);
  //if (rc != 0) ERR_MSG("setvbuf %s\n",strerror(errno));

  rc = sem_init(&(ctx->task),0,0);
  if (rc != 0) ERR_MSG("sem_init %s\n",strerror(errno));

  ctx->num_thread_groups = num_thread_groups;

  ctx->data = calloc(ctx->num_thread_groups,sizeof(float*));

  ctx->dataoffset = calloc(ctx->num_thread_groups,sizeof(size_t));
  ctx->datasize   = calloc(ctx->num_thread_groups,sizeof(size_t));

  size_t fwdsize = 1LL * nnx * nny * diam_width;

  for (int i = 0; i < ctx->num_thread_groups; i ++) {
    ctx->data[i] = calloc(fwdsize,sizeof(float));
  }

  ctx->index = calloc(ctx->num_thread_groups,sizeof(int));
  if (ctx->index == NULL) ERR_MSG("calloc %s\n",strerror(errno));
  ctx->inuse = calloc(ctx->num_thread_groups,sizeof(int));
  if (ctx->inuse == NULL) ERR_MSG("calloc %s\n",strerror(errno));

  ctx->head = 0;
  ctx->tail = 0;

  ctx->mutex_p = calloc(ctx->num_thread_groups,sizeof(pthread_mutex_t));
  ctx->mutex_c = calloc(ctx->num_thread_groups,sizeof(pthread_mutex_t));
  ctx->mutex_cond = calloc(ctx->num_thread_groups,sizeof(pthread_cond_t));
  rc = pthread_mutex_init(&(ctx->mutex),NULL);
  if (rc != 0) ERR_MSG("pthread_mutex_init %s\n",strerror(errno));

  for (int i = 0; i < ctx->num_thread_groups; i ++) {
     rc = pthread_mutex_init(&(ctx->mutex_p[i]),NULL);
     if (rc != 0) ERR_MSG("pthread_mutex_init %s\n",strerror(errno));

     rc = pthread_mutex_init(&(ctx->mutex_c[i]),NULL);
     if (rc != 0) ERR_MSG("pthread_mutex_init %s\n",strerror(errno));

     rc = pthread_cond_init(&(ctx->mutex_cond[i]),NULL);
     if (rc != 0) ERR_MSG("pthread_cond_init %s\n",strerror(errno));
  }

  ctx->coreid = coreid;

  rc = pthread_create(&(ctx->thread_id), NULL, mlbs_writer_helper, ctx);
  if (rc != 0) ERR_MSG("pthread_create %s\n",strerror(errno));

  return 0;
}

int mlbs_write(mlbs_t * ctx,
               int groupid,
               float  * data,
               size_t dataoffset,
               size_t datasize) {
  int rc;

  // check groupid is free
  if (groupid != -1) {
    rc = pthread_mutex_lock(&ctx->mutex_c[groupid]);
    if (rc != 0) ERR_MSG("pthread_mutex_lock %s\n",strerror(errno));

    if (ctx->inuse[groupid] == 1) {
      rc = pthread_cond_wait(&ctx->mutex_cond[groupid],&ctx->mutex_c[groupid]);
      if (rc != 0) ERR_MSG("pthread_cond_wait %s\n",strerror(errno));
    }

    rc = pthread_mutex_unlock(&ctx->mutex_c[groupid]);
    if (rc != 0) ERR_MSG("pthread_mutex_lock %s\n",strerror(errno));

    memcpy(ctx->data[groupid],data,datasize*sizeof(float));
  }
  // enqueue data write
  rc = pthread_mutex_lock(&ctx->mutex);
  if (rc != 0) ERR_MSG("pthread_mutex_lock %s\n",strerror(errno));

  if (groupid != -1) {
    ctx->inuse[groupid] = 1;
  }
  ctx->index[ctx->tail]      = groupid;
  ctx->datasize[ctx->tail]   = datasize;
  ctx->dataoffset[ctx->tail] = dataoffset;

  ctx->tail ++;
  if (ctx->tail == ctx->num_thread_groups) ctx->tail = 0;

  rc = pthread_mutex_unlock(&ctx->mutex);
  if (rc != 0) ERR_MSG("pthread_mutex_unlock %s\n",strerror(errno));

  rc  = sem_post(&ctx->task);
  if (rc != 0) ERR_MSG("sem_post %s\n",strerror(errno));
}

void * mlbs_writer_helper(void * data) {
  int rc, groupid;
  size_t offset;
  mlbs_t * ctx = (mlbs_t *) data;
  size_t size_done,datasize;

  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(ctx->coreid, &cpuset);
  rc = sched_setaffinity(0,sizeof(cpu_set_t), &cpuset);

  while (1) {

    rc = sem_wait(&ctx->task);
    if (rc != 0) ERR_MSG("sem_wait %s\n",strerror(errno));

    rc = pthread_mutex_lock(&ctx->mutex);
    if (rc != 0) ERR_MSG("pthread_mutex_lock %s\n",strerror(errno));

    groupid = ctx->index[ctx->head];
    if (groupid == -1) {
      rc = pthread_mutex_unlock(&ctx->mutex);
      if (rc != 0) ERR_MSG("pthread_mutex_unlock %s\n",strerror(errno));
      break;
    }

    offset   = ctx->dataoffset[ctx->head];
    datasize = ctx->datasize[ctx->head];

    ctx->head ++;
    if (ctx->head == ctx->num_thread_groups) ctx->head = 0;


    rc = pthread_mutex_unlock(&ctx->mutex);
    if (rc != 0) ERR_MSG("pthread_mutex_unlock %s\n",strerror(errno));

    rc = fseek(ctx->fp, offset, SEEK_SET);
    if (rc != 0) ERR_MSG("fseek %s\n",strerror(errno));

    size_done = fwrite(ctx->data[groupid],sizeof(float),datasize,ctx->fp);
    if (size_done != datasize) ERR_MSG("fwrite %s %d %d\n",strerror(errno),size_done,datasize);

    // signal the slot is available

    rc = pthread_mutex_lock(&ctx->mutex_c[groupid]);
    if (rc != 0) ERR_MSG("pthread_mutex_lock %s\n",strerror(errno));

    ctx->inuse[groupid] = 0;

    rc = pthread_cond_signal(&ctx->mutex_cond[groupid]);
    if (rc != 0) ERR_MSG("pthread_cond_signal %s\n",strerror(errno));

    rc = pthread_mutex_unlock(&ctx->mutex_c[groupid]);
    if (rc != 0) ERR_MSG("pthread_mutex_lock %s\n",strerror(errno));

  }

  printf("END OF WRITER HELPER\n");
}

int mlbs_free(mlbs_t * ctx) {

  int rc;

  rc = pthread_join(ctx->thread_id, NULL);
  if (rc != 0) ERR_MSG("pthread_join %s\n",strerror(errno));

  rc = sem_destroy(&(ctx->task));
  if (rc != 0) ERR_MSG("sem_destroy %s\n",strerror(errno));

  for (int i = 0; i < ctx->num_thread_groups; i ++) {
    free(ctx->data[i]); ctx->data[i] = NULL;
  }

  free(ctx->data); ctx->data = NULL;

  free(ctx->index); ctx->index = NULL;
  free(ctx->inuse); ctx->inuse = NULL;

  rc = pthread_mutex_destroy(&(ctx->mutex));
  if (rc != 0) ERR_MSG("pthread_mutex_destroy %s\n",strerror(errno));

  for (int i = 0; i < ctx->num_thread_groups; i ++) {
    rc = pthread_mutex_destroy(&(ctx->mutex_p[i]));
    if (rc != 0) ERR_MSG("pthread_mutex_destroy %s\n",strerror(errno));

    rc = pthread_mutex_destroy(&(ctx->mutex_c[i]));
    if (rc != 0) ERR_MSG("pthread_mutex_destroy %s\n",strerror(errno));

    rc = pthread_cond_destroy(&(ctx->mutex_cond[i]));
    if (rc != 0) ERR_MSG("pthread_cond_destroy %s\n",strerror(errno));

  }

  free(ctx->mutex_p); ctx->mutex_p = NULL;
  free(ctx->mutex_c); ctx->mutex_c = NULL;
  free(ctx->mutex_cond); ctx->mutex_cond = NULL;

  free(ctx->dataoffset); ctx->dataoffset = NULL;
  free(ctx->datasize);   ctx->datasize   = NULL;

  fclose(ctx->fp);
  return 0;
}