///
/// @copyright Copyright . All rights reserved.
/// This file is part of \b stencil.
///
/// @b stencil is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// @b stencil is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with \b stencil.  If not, see <http://www.gnu.org/licenses/>.
///
/// @author
/// @file src/gpu_wave.cu
/// @brief This file contains the implementation of the new GPU wave descriptor.
///
#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h>
#include <math.h>
#include <stencil/config.h>
#include <stencil/macros.h>
#include <stencil/sismap.h>
#include <stencil/shot.h>

/// @brief This function should contains the OpenCL
/// kernel that updates the wave fields on the GPU.
/// Note: do not change the signature of the function
/// unless you use local memory 
/// @param u0 the wave fields at time step \em t
/// @param u1 the wave fields at time step \em t+1
/// @param roc2 the velocity terms
/// @param eta PML initial terms
/// @param phi PML altered terms
/// @param coefx the stencil coefficients on the \b x axis
/// @param coefy the stencil coefficients on the \b y axis
/// @param coefz the stencil coefficients on the \b z axis
/// @param dimx the domain width
/// @param dimy the domain height
/// @param dimz the domain depth
/// @param sx the stencil half length on the \b x axis
/// @param sy the stencil half length on the \b y axis
/// @param sz the stencil half length on the \b z axis
/// @param pmlx the PML width on the \b x axis
/// @param pmly the PML width on the \b y axis
/// @param pmlz the PML width on the \b z axis
/// @param hdx2 the PML absorbing factor on the \b x axis
/// @param hdy2 the PML absorbing factor on the \b y axis
/// @param hdz2 the PML absorbing factor on the \b z axis
#define   U0(x,y,z) (  u0[(2*sx+dimx)*((2*sy+dimy)*(z+sz) + y+sy) + x+sx])
#define   U1(x,y,z) (  u1[(2*sx+dimx)*((2*sy+dimy)*(z+sz) + y+sy) + x+sx])
#define ROC2(x,y,z) (roc2[dimx*(dimy*(z) + y ) + x]) 
#define  ETA(x,y,z) ( eta[(dimx+2)*((dimy+2)*(z+1) + (y+1) ) + (x+1)]) 
#define  PHI(x,y,z) (phi[dimx*(dimy*(z) + y ) + x])
#define  IMGT(tab, x,y,z) (tab[x + (img_dimx * (img_dimy*(z) + y))])
#define  FWD(x,y,z) ( fwd[(2*sx+dimx)*((2*sy+dimy)*(z+sz) + y+sy) + x+sx])

__global__  void cu_wave_update_fields_inner(float          *u0,   
                                             float          *u1,
                                             float        *roc2,
                                             const float *coefx,
                                             const float *coefy,
                                             const float *coefz,
                                             unsigned int  dimx, 
                                             unsigned int  dimy, 
                                             unsigned int  dimz,
                                             int             sx,   
                                             int             sy, 
                                             int             sz,
                                             int xbeg, int xend,
                                             int ybeg, int yend,
                                             int zbeg, int zend) {
  /// global idx
  int xgid = threadIdx.x + blockIdx.x*blockDim.x;
  int ygid = threadIdx.y + blockIdx.y*blockDim.y;

  int zgid;

  /// laplacian operator
  float laplacian;
  float   current; 
  float  front[4];
  float behind[4];
  float coef0 = coefx[0] + coefy[0] + coefz[0];

  if( (xgid>=xbeg) && (xgid<xend) && (ygid>=ybeg) && (ygid<yend) ) { 
    current   = 0.0;  
    behind[0] = 0.0;
    behind[1] = 0.0;
    behind[2] = 0.0;
    behind[3] = 0.0;

    /// for register blocking 
    front[0] = U0(xgid,ygid,0);
    front[1] = U0(xgid,ygid,1);
    front[2] = U0(xgid,ygid,2);
    front[3] = U0(xgid,ygid,3);

    for(zgid = zbeg; zgid < zend; zgid++) {
      behind[3] = behind[2];                   
      behind[2] = behind[1];                   
      behind[1] = behind[0];                   
      behind[0] = current;                    
      current   = front[0];                    
      front[0]  = front[1];                    
      front[1]  = front[2];                    
      front[2]  = front[3];                    
      front[3]  = U0(xgid, ygid, (zgid+sz));

      laplacian  = coef0 * current;                         
      laplacian += coefx[1] * ( U0(xgid+1, ygid, zgid) +    
        U0(xgid-1, ygid, zgid) );   
      laplacian += coefx[2] * ( U0(xgid+2, ygid, zgid) +    
        U0(xgid-2, ygid, zgid) );   
      laplacian += coefx[3] * ( U0(xgid+3, ygid, zgid) +    
        U0(xgid-3, ygid, zgid) );   
      laplacian += coefx[4] * ( U0(xgid+4, ygid, zgid) +    
        U0(xgid-4, ygid, zgid) );   

      laplacian += coefy[1] * ( U0(xgid, ygid+1, zgid) +    
        U0(xgid, ygid-1, zgid) );   
      laplacian += coefy[2] * ( U0(xgid, ygid+2, zgid) +    
        U0(xgid, ygid-2, zgid) );   
      laplacian += coefy[3] * ( U0(xgid, ygid+3, zgid) +    
        U0(xgid, ygid-3, zgid) );   
      laplacian += coefy[4] * ( U0(xgid, ygid+4, zgid) +    
        U0(xgid, ygid-4, zgid) );   

      laplacian += coefz[1] * ( front[0] + behind[0] );       
      laplacian += coefz[2] * ( front[1] + behind[1] );       
      laplacian += coefz[3] * ( front[2] + behind[2] );       
      laplacian += coefz[4] * ( front[3] + behind[3] );       

      U1(xgid, ygid, zgid) = 2.0 * current - U1(xgid, ygid, zgid) + 
                             ROC2(xgid, ygid, zgid)*laplacian;
    }
  }
}

__global__  void cu_wave_update_fields(float          *u0,   
                                       float          *u1,
                                       float        *roc2,
                                       float         *eta,
                                       float         *phi,
                                       const float *coefx,
                                       const float *coefy,
                                       const float *coefz,
                                       unsigned int  dimx, 
                                       unsigned int  dimy, 
                                       unsigned int  dimz,
                                       int             sx,   
                                       int             sy, 
                                       int             sz,
                                       int           pmlx,   
                                       int           pmly, 
                                       int           pmlz,
                                       float         hdx2,
                                       float         hdy2, 
                                       float         hdz2) {
  /// global idx
  int xgid = threadIdx.x + blockIdx.x*blockDim.x;
  int ygid = threadIdx.y + blockIdx.y*blockDim.y;

  int zgid;

  /// laplacian operator
  float laplacian;
  float   current; 
  float  front[4];
  float behind[4];
  float coef0 = coefx[0] + coefy[0] + coefz[0];

  if( (xgid<dimx) && (ygid<dimy) ) { 
    current   = 0.0;  
    behind[0] = 0.0;
    behind[1] = 0.0;
    behind[2] = 0.0;
    behind[3] = 0.0;

    /// for register blocking 
    front[0] = U0(xgid,ygid,0);
    front[1] = U0(xgid,ygid,1);
    front[2] = U0(xgid,ygid,2);
    front[3] = U0(xgid,ygid,3);

    for(zgid = 0; zgid < dimz; zgid++) {
      behind[3] = behind[2];                   
      behind[2] = behind[1];                   
      behind[1] = behind[0];                   
      behind[0] = current;                    
      current   = front[0];                    
      front[0]  = front[1];                    
      front[1]  = front[2];                    
      front[2]  = front[3];                    
      front[3]  = U0(xgid, ygid, (zgid+sz));

      laplacian  = coef0 * current;                         
      laplacian += coefx[1] * ( U0(xgid+1, ygid, zgid) +    
        U0(xgid-1, ygid, zgid) );   
      laplacian += coefx[2] * ( U0(xgid+2, ygid, zgid) +    
        U0(xgid-2, ygid, zgid) );   
      laplacian += coefx[3] * ( U0(xgid+3, ygid, zgid) +    
        U0(xgid-3, ygid, zgid) );   
      laplacian += coefx[4] * ( U0(xgid+4, ygid, zgid) +    
        U0(xgid-4, ygid, zgid) );   

      laplacian += coefy[1] * ( U0(xgid, ygid+1, zgid) +    
        U0(xgid, ygid-1, zgid) );   
      laplacian += coefy[2] * ( U0(xgid, ygid+2, zgid) +    
        U0(xgid, ygid-2, zgid) );   
      laplacian += coefy[3] * ( U0(xgid, ygid+3, zgid) +    
        U0(xgid, ygid-3, zgid) );   
      laplacian += coefy[4] * ( U0(xgid, ygid+4, zgid) +    
        U0(xgid, ygid-4, zgid) );   

      laplacian += coefz[1] * ( front[0] + behind[0] );       
      laplacian += coefz[2] * ( front[1] + behind[1] );       
      laplacian += coefz[3] * ( front[2] + behind[2] );       
      laplacian += coefz[4] * ( front[3] + behind[3] );       

      U1(xgid, ygid, zgid) =
      (
       (2.-ETA(xgid, ygid, zgid)*ETA(xgid, ygid, zgid) + 2.*ETA(xgid, ygid, zgid))
       * current
       - U1(xgid, ygid, zgid)
       + ROC2(xgid, ygid, zgid) * (laplacian + PHI(xgid, ygid, zgid))
       ) / (1. + 2.*ETA(xgid, ygid, zgid)); 

      PHI(xgid, ygid, zgid) =                 
      (
       PHI(xgid, ygid, zgid)-                                     
       ((ETA(xgid+1, ygid,     zgid) - ETA(xgid-1, ygid,     zgid))      
         *( U0(xgid+1, ygid,     zgid) -  U0(xgid-1, ygid,     zgid))*hdx2 
         +(ETA(xgid,   ygid+1,   zgid) - ETA(xgid,   ygid-1,   zgid))      
         *( U0(xgid,   ygid+1,   zgid) -  U0(xgid,   ygid-1,   zgid))*hdy2 
         +(ETA(xgid,   ygid,   zgid+1) - ETA(xgid,   ygid,   zgid-1))      
         *( U0(xgid,   ygid,   zgid+1) -  U0(xgid,   ygid,   zgid-1))*hdz2)
       ) / (1. + ETA(xgid, ygid, zgid));
    }
  }
}

/// @brief This function should contains the OpenCL
/// kernel that updates the wave source on the GPU.
/// Note: this kernel can be omitted
/// 
/// @param u1 the wave fields at time step \em t+1
/// @param src the array that contains the source terms
/// @param dimx the domain width
/// @param dimy the domain height
/// @param sx the stencil half length on the \b x axis
/// @param sy the stencil half length on the \b y axis
/// @param sz the stencil half length on the \b z axis
/// @param srcx the source location on the \b x axis
/// @param srcy the source location on the \b y axis
/// @param srcz the source location on the \b z axis
/// @param t the current time step
__global__  void cu_wave_update_source(float*         u1,
                                       float       sterm,
                                       unsigned int dimx,
                                       unsigned int dimy, 
                                       int            sx,
                                       int            sy,
                                       int            sz,
                                       int        srcidx, 
                                       int         depth) {
  u1[(2*sx+dimx)*(2*sy+dimy)*(depth+sz) + srcidx] += sterm;
}

__global__  void cu_wave_extract_sismos(float          *u1,
                                        float      *sismos,
                                        unsigned int    *r,
                                        unsigned int  dimx, 
                                        unsigned int  dimy, 
                                        int             sx,
                                        int             sy,
                                        int             sz,
                                        unsigned int     n,
                                        unsigned int depth, 
                                        unsigned int t) {
  /// global index.
  int ir = threadIdx.x + blockIdx.x*blockDim.x;
  if (ir < n)
    sismos[ir + n*t] = u1[(2*sx+dimx)*(2*sy+dimy)*(depth+sz) + r[ir]];
}

__global__  void cu_wave_inject_sismos(float          *u0,
                                       float      *sismos,
                                       unsigned int    *r,
                                       unsigned int  dimx, 
                                       unsigned int  dimy, 
                                       int             sx,
                                       int             sy,
                                       int             sz,
                                       unsigned int     n,
                                       unsigned int depth, 
                                       unsigned int t) {
  /// global index.
  int ir = threadIdx.x + blockIdx.x*blockDim.x;
  if (ir < n)
    u0[(2*sx+dimx)*(2*sy+dimy)*(depth+sz) + r[ir]] += sismos[ir + n*t];
}

__global__  void cu_wave_img_cond(float          *u1,   
                                  float         *fwd,
                                  float         *img,
                                  float         *ilm,
                                  unsigned int  img_dimx, 
                                  unsigned int  img_dimy, 
                                  unsigned int  img_dimz,
                                  unsigned int  dimx, 
                                  unsigned int  dimy, 
                                  unsigned int  dimz,
                                  int             sx,   
                                  int             sy, 
                                  int             sz,
                                  int px, int py) {
  /// global index.
  int xgid = threadIdx.x + blockIdx.x*blockDim.x;
  int ygid = threadIdx.y + blockIdx.y*blockDim.y;
  if( (xgid<img_dimx) && (ygid<img_dimy) ) { 
    for(int zgid = 0; zgid < img_dimz; zgid++) {
      IMGT(img, xgid, ygid, zgid) +=  
        U1(xgid+px,ygid+py,zgid)*FWD(xgid+px,ygid+py,zgid);
      IMGT(ilm, xgid, ygid, zgid) += 
        FWD(xgid+px,ygid+py,zgid)*FWD(xgid+px,ygid+py,zgid);
    }
  }
}

__global__  void cu_wave_img_gather(float         *img,
                                    float         *img_shot,
                                    float         *ilm_shot,
                                    unsigned int  img_dimx, 
                                    unsigned int  img_dimy, 
                                    unsigned int  img_dimz,
                                    unsigned int  nb_shots) {
  /// global index.
  int xgid = threadIdx.x + blockIdx.x*blockDim.x;
  int ygid = threadIdx.y + blockIdx.y*blockDim.y;
  if( (xgid<img_dimx) && (ygid<img_dimy) ) { 
    for(int zgid = 0; zgid < img_dimz; zgid++) {
      IMGT(img, xgid, ygid, zgid) += IMGT(img_shot, xgid, ygid, zgid);
        ///(IMGT(ilm_shot, xgid, ygid, zgid)); 
    }
  }
}

extern "C" void gpu_wave_set(int device) {
  int nb_devices;
  GPU_CHK(cudaGetDeviceCount(&nb_devices));
  ERR_IF(nb_devices == 0, "cannot find any GPU device");
  ERR_IF(device >= nb_devices, "invalid device index");
  cudaSetDevice(device);
}

extern "C" void gpu_wave_unset() {
  cudaDeviceReset();
}

extern "C" void gpu_wave_init(sismap_t *s) {
  GPU_CHK(cudaMalloc((void**)&s->d_coefx, (s->sx+1)*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&s->d_coefy, (s->sy+1)*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&s->d_coefz, (s->sz+1)*sizeof(float)));
  GPU_CHK(cudaMalloc((void**)&s->d_rcv, (s->rcv_len)*sizeof(unsigned int)));

  GPU_CHK(cudaMemcpy(s->d_coefx, s->coefx, 
                     (s->sx+1)*sizeof(float),
                     cudaMemcpyHostToDevice));
  GPU_CHK(cudaMemcpy(s->d_coefy, s->coefy, 
                     (s->sy+1)*sizeof(float),
                     cudaMemcpyHostToDevice));
  GPU_CHK(cudaMemcpy(s->d_coefz, s->coefz, 
                     (s->sz+1)*sizeof(float),
                     cudaMemcpyHostToDevice));
  GPU_CHK(cudaMemcpy(s->d_rcv, s->rcv, 
                    (s->rcv_len)*sizeof(unsigned int),
                     cudaMemcpyHostToDevice));
  s->l.x = 32;
  s->l.y = 8;
  s->l.z = 1;
  s->g.x = (s->dimx + s->l.x -1)/s->l.x;
  s->g.y = (s->dimy + s->l.y -1)/s->l.y;
  s->g.z = 1;
  s->g_img.x = (s->img_dimx + s->l.x -1)/s->l.x;
  s->g_img.y = (s->img_dimy + s->l.y -1)/s->l.y;
  s->g_img.z = 1;
  s->o.x = s->o.y = s->o.z = 1;
}

extern "C" void gpu_wave_release(sismap_t *s) {
  GPU_CHK(cudaFree(s->d_coefx));
  GPU_CHK(cudaFree(s->d_coefy));
  GPU_CHK(cudaFree(s->d_coefz));
  GPU_CHK(cudaFree(s->d_rcv));
}

extern "C" void gpu_wave_info(sismap_t *s) {
  cudaDeviceProp props;
  /// query device properties.
  GPU_CHK(cudaGetDeviceProperties(&props, s->device));
  MSG("... GPU device: %s", props.name);
}

extern "C" void gpu_wave_update_fields(sismap_t *s, 
                                      float* d_u0, float *d_u1, 
                                      float *d_roc2,
                                      float* phi, float* eta) {
  cu_wave_update_fields<<<s->g, s->l>>>(d_u0, d_u1, d_roc2,eta,
                                                           phi,
                                                           s->d_coefx,
                                                           s->d_coefy,
                                                           s->d_coefz,
                                                           s->dimx,
                                                           s->dimy,
                                                           s->dimz,
                                                           s->sx,
                                                           s->sy,
                                                           s->sz,
                                                           s->pmlx,
                                                           s->pmly,
                                                           s->pmlz,
                                                           s->hdx2,
                                                           s->hdy2,
                                                           s->hdz2);
  GPU_CHK(cudaGetLastError());
}

extern "C" void gpu_wave_update_source(sismap_t *s, 
                                       shot_t *shot, float *d_u0, float sterm) {
  cu_wave_update_source<<<s->o, s->o>>>(d_u0, sterm,
                                        s->dimx,
                                        s->dimy,
                                        s->sx,
                                        s->sy,
                                        s->sz,
                                        shot->srcidx,
                                        s->src_depth);
  GPU_CHK(cudaGetLastError());
}

extern "C" void gpu_wave_extract_sismos(sismap_t *s, 
                                        float* d_u0, 
                                        unsigned int t, float *d_sismos) {
  cu_wave_extract_sismos<<<(s->rcv_len + 127)/128, 128>>>(d_u0, d_sismos,
                                                           s->d_rcv,
                                                           s->dimx,
                                                           s->dimy,
                                                           s->sx,
                                                           s->sy,
                                                           s->sz,
                                                           s->rcv_len,
                                                           s->rcv_depth, t);
  GPU_CHK(cudaGetLastError());
}

extern "C" void gpu_wave_inject_sismos(sismap_t *s, 
                                       float* d_u0, 
                                       unsigned int t, float *d_sismos) {
  cu_wave_inject_sismos<<<(s->rcv_len + 127)/128, 128>>>(d_u0, d_sismos,
                                                           s->d_rcv,
                                                           s->dimx,
                                                           s->dimy,
                                                           s->sx,
                                                           s->sy,
                                                           s->sz,
                                                           s->rcv_len,
                                                           s->rcv_depth, t);
  GPU_CHK(cudaGetLastError());
}

extern "C" void gpu_wave_image_condition(sismap_t* s, 
                                         float *d_u1, float *d_fwd, 
                                         float *d_img, 
                                         float *d_ilm, unsigned int t) {
  if(t%s->nb_snap==0) {
    cu_wave_img_cond<<<s->g_img, s->l>>>(d_u1,
                                         d_fwd,
                                         d_img,
                                         d_ilm,
                                         s->img_dimx,
                                         s->img_dimy,
                                         s->img_dimz,
                                         s->dimx,
                                         s->dimy,
                                         s->dimz,
                                         s->sx,
                                         s->sy,
                                         s->sz, s->pmlx, s->pmly);
    GPU_CHK(cudaGetLastError());
    //GPU_CHK(cudaDeviceSynchronize());
  }
}

extern "C" void gpu_wave_image_gather(sismap_t* s, 
                                      float *d_img, 
                                      float *d_ilm, float *gimg) {
  cu_wave_img_gather<<<s->g_img, s->l>>>(gimg,
                                         d_img,
                                         d_ilm,
                                         s->img_dimx,
                                         s->img_dimy,
                                         s->img_dimz, s->nb_shots);
  GPU_CHK(cudaGetLastError());
}
/// @brief An access helper to s->u0
#define HU0(z,y,x)                                                    \
(u0[(x+s->sx) + (2*s->sx + s->dimx)*((2*s->sy + s->dimy) * (z+s->sz)  \
   + (y+s->sy))])

/// @brief An access helper to s->u1
#define HU1(z,y,x)                   \
(u1[(x+s->sx) + (2*s->sx + s->dimx) *\
  ((2*s->sy + s->dimy) * (z+s->sz) + (y+s->sy))])

/// @brief Note that the matrices are transposed.
/// This is for visualization purposes only.
/// Please do not change this behavior. However, it is 
/// recommended to overlay communication with computation.
extern "C" void gpu_wave_save_snapshot(sismap_t *s, shot_t *shot,
                                       float *d_u1, float *u1, unsigned int t) {
  if (t%s->nb_snap==0) {
    GPU_CHK(cudaMemcpy(u1, d_u1, 
                       s->size*sizeof(float),
                       cudaMemcpyDeviceToHost));
    CHK(fwrite(u1, sizeof(float), s->size, shot->fd_snap) != s->size,
        "failed to write snapshot");
    if (s->verbose) MSG("... GPU saving snapshot %4d t=%u", s->snap_idx, t);
    s->snap_idx=s->snap_idx+1;
  }
}

extern "C" void gpu_wave_read_snapshot(sismap_t* s, 
                                       shot_t *shot, 
                                       float *d_fwd, float *fwd, unsigned int t) {
  if(t%s->nb_snap==0) {
    if (s->verbose) MSG("... GPU reading snapshot %4d t=%u", s->snap_idx, t);    
    CHK(fseek(shot->fd_snap, s->snap_idx*s->size*sizeof(float), SEEK_SET) != 0,
        "failed to fseek file");
    CHK(fread(fwd, sizeof(float), s->size, shot->fd_snap) != s->size, 
        "failed to read snapshot");
    s->snap_idx = s->snap_idx-1;
    GPU_CHK(cudaMemcpy(d_fwd, fwd, 
                       s->size*sizeof(float), cudaMemcpyHostToDevice));
    //GPU_CHK(cudaDeviceSynchronize());
  }
}

extern "C" unsigned int gpu_wave_check_fields(sismap_t *s, 
                                              float *d_u0, float *u1, 
                                              float *u0,  float e) {
  unsigned int x, y, z, not_valid = 0;
  memset(u1, 0, s->size*sizeof(float));
  GPU_CHK(cudaMemcpy(u1, d_u0, 
                     s->size*sizeof(float),
                     cudaMemcpyDeviceToHost));
  for(z = 0; z < s->dimz; ++z) {
    for(y = 0; y < s->dimy; ++y) {
      for(x = 0; x < s->dimx; ++x) {
        if(isnan(HU1(z,y,x)) || (fabs(HU0(z,y,x) - HU1(z,y,x))>e)) {
          not_valid++;
          #ifdef __DEBUG
          fprintf(stderr, "%3d %3d %3d -> cpu: %-12.6f, gpu: %-12.6f\n",
                  z, y, x, HU0(z,y,x), HU1(z,y,x));
          #endif // __DEBUG
        }
      }
    }
  }
  return not_valid;
}

extern "C" void gpu_wave_save_sismos(sismap_t* s, 
                                     shot_t *shot, 
                                     float* d_sismos, float *sismos) {
  size_t sz = s->rcv_len*s->time_steps*sizeof(float); 
  GPU_CHK(cudaMemcpy(sismos, d_sismos, sz, cudaMemcpyDeviceToHost));
  char tmp[512];
  sprintf(tmp, "%s/gpu_%s_%d.raw", OUTDIR, SISMOS_BASE, shot->id);
  FILE *fd = fopen(tmp, "wb");
  CHK(fd == NULL, "failed to open sismos file");
  CHK(fwrite(sismos, sz, 1, fd) != 1, "failed to properly write sismos");
  fclose(fd);
}

extern "C" void gpu_wave_read_sismos(sismap_t* s, 
                                     shot_t *shot, 
                                     float* d_sismos, float *sismos) {
  char tmp[512];
  FILE *fd;
  size_t sz = s->rcv_len*s->time_steps*sizeof(float); 
  sprintf(tmp, "%s/gpu_%s_%d.raw", OUTDIR, SISMOS_BASE, shot->id);
  fd = fopen(tmp, "rb");
  CHK(fd == NULL, 
      "sismos file not found, aborting (run modeling to generate it)"); 
  CHK(fread(sismos, sz, 1, fd) != 1, "failed to properly read sismos");
  fclose(fd);
  GPU_CHK(cudaMemcpy(d_sismos, sismos, sz, cudaMemcpyHostToDevice));
}

extern "C" void gpu_wave_save_fwd_dbg(sismap_t* s, 
                                      shot_t *shot, 
                                      float *d_u1, float* u1, unsigned int t) {
  #ifdef __DEBUG
  if(t%s->nb_snap==0) {
    GPU_CHK(cudaMemcpy(u1, d_u1, 
                       s->size*sizeof(float),
                       cudaMemcpyDeviceToHost));
    GPU_CHK(cudaDeviceSynchronize());
    for (unsigned int z = 0; z < s->dimz; z++)
      for (unsigned int x = 0; x < s->dimx; x++)
        fwrite(&HU1(z, shot->srcidx/(s->dimx+2*s->sx), x), 
               sizeof(float), 1, shot->fd_fwd);
  }
  #endif // __DEBUG
}

extern "C" void gpu_wave_save_bwd_dbg(sismap_t* s, 
                                      shot_t *shot, 
                                      float *d_u1, float* u1, unsigned int t) {
  #ifdef __DEBUG
  if(t%s->nb_snap==0) {
    GPU_CHK(cudaMemcpy(u1, d_u1, 
                       s->size*sizeof(float),
                       cudaMemcpyDeviceToHost));
    GPU_CHK(cudaDeviceSynchronize());
    for (unsigned int z = 0; z < s->dimz; z++)
      for (unsigned int x = 0; x < s->dimx; x++)
        fwrite(&HU1(z, shot->srcidx/(s->dimx+2*s->sx), x), 
                sizeof(float), 1, shot->fd_bwd);
  }
  #endif // __DEBUG
}

#define HIMG(z,y,x)  (img[x + (s->img_dimx * (s->img_dimy*(z) + y))])

extern "C" void gpu_wave_save_img(sismap_t* s, shot_t *shot, 
                                  float *d_img, float *d_ilm, float* tmp) {
  // cu_wave_img_gather<<<s->g_img, s->l>>>(d_img,
  //                                        d_img,
  //                                        d_ilm,
  //                                        s->img_dimx,
  //                                        s->img_dimy,
  //                                        s->img_dimz, s->nb_shots);
  // GPU_CHK(cudaGetLastError());
  GPU_CHK(cudaMemcpy(tmp, d_ilm, 
                     s->size_img*sizeof(float),
                     cudaMemcpyDeviceToHost));
  GPU_CHK(cudaDeviceSynchronize());
  fwrite(tmp, sizeof(float), s->size_img, shot->fd_ilm);
  GPU_CHK(cudaMemcpy(tmp, d_img, 
                     s->size_img*sizeof(float),
                     cudaMemcpyDeviceToHost));
  GPU_CHK(cudaDeviceSynchronize());
  fwrite(tmp, sizeof(float), s->size_img, shot->fd_img);
}

/// @brief undefine static macros
#undef HU0
#undef HU1
#undef GPU_CHK
#undef HIMG
