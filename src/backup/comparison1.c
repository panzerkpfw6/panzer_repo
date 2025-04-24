/*

void femwd_iso_ref_1st_orig( const int shape[3], const int zb, const int yb_r0,
                        const int xb, const int ze, const int ye_r0, const int xe,
                    const real_t *  coef, hFloat *  p11, hFloat *  p12, hFloat *  p13,
                    hFloat *  p21, hFloat *  p22, hFloat *  p23,const hFloat *  roc2,
                    float * dampx,float * dampy,float * dampz,
                    int t_dim, int b_inc, int e_inc,int NHALO,
                    int tb, int te,int t0, stencil_ctx stencil_ctx,int mtid,tb_data_t * data)
{
#pragma omp parallel shared(shape, stencil_ctx, roc2, coef, mtid, tb, te, t_dim, NHALO,recv_rec,irecv_rec) \
firstprivate(b_inc, e_inc) \
num_threads(stencil_ctx.thread_group_size)
    {
        int lstencil=NHALO;// @pavel  allocate variable lstencil
        int tgs, nwf, th_nwf, tid, gtid, xi, yb, ye, ib, ie, kt, t,  q, r, err;
        double t_start;

        const int nnx =shape[2];
        const int nny =shape[1];
        const int nnz =shape[0];
        const unsigned long nnzy = 1UL * nnz * nny;
        const unsigned long nnyz = nnzy;
        // index notation for velocity array
        const int nnz_v=stencil_ctx.nz;
        const unsigned long nnyz_v=1UL*stencil_ctx.nz*stencil_ctx.ny;

        tgs = stencil_ctx.thread_group_size;
        nwf = stencil_ctx.num_wf;

        tid = 0;
        gtid = 0;
#if defined(_OPENMP)
        tid = omp_get_thread_num();
		gtid = tid + mtid * tgs;
#endif


        if(stencil_ctx.use_manual_cpu_bind == 1){
            err = sched_setaffinity(0, stencil_ctx.setsize, stencil_ctx.bind_masks[mtid*tgs+tid]);
            if(err==-1) printf("WARNING: Could not set CPU Affinity\n");
        }

        hFloat *  u1 = p11;
        hFloat *  u2 = p12;
        hFloat *  u3 = p13;
        hFloat *  v1 = p21;
        hFloat *  v2 = p22;
        hFloat *  v3 = p23;


        int th_z = stencil_ctx.th_z;
        int th_y = stencil_ctx.th_y;
        int th_x = stencil_ctx.th_x;

        // tid = tid_x*(th_z*th_y) + tid_y*th_z + tid_z
        int tid_z = tid%th_z;
        int tid_y = tid/th_z;
        int tid_x = tid/(th_z*th_y);

        int yb_r = yb_r0;
        int ye_r = ye_r0;

        if(stencil_ctx.th_y>1 ){
            if(b_inc !=0 && e_inc!=0){ // split only at full diamonds
                if (tid_y%2 == 0){ // left thread
                    ye_r = (yb_r + ye_r)/2;
                    e_inc = 0;
                } else{
                    yb_r = (yb_r + ye_r)/2;
                    b_inc = 0;
                }
            }else{// use the y-threads along x-axis make sure to use sufficient number of frontlines
                th_x *= th_y;
                tid_x = tid/th_z;
                if (nwf < th_x) nwf = th_x;
            }
        }

        int nbz = (ze-zb)/th_z;
        q = (int)((ze-zb)/th_z);
        r = (ze-zb)%th_z;
        if(tid_z < r) {
            ib = zb + tid_z * (q+1);
            ie = ib + (q+1);
        }else {
            ib = zb + r * (q+1) + (tid_z - r) * q;
            ie = ib + q;
        }

        th_nwf = nwf/th_x;

        int printed = 0; //@KADIR
        int iz_=data->rcv_depth; //@pavel
        int end=0;

        const Myfloat inv_dx = 1. / (stencil_ctx.dx);
        const Myfloat inv_dy = 1. / (stencil_ctx.dy);
        const Myfloat inv_dz = 1. / (stencil_ctx.dz);

        for(xi=xb; xi<xe; xi+=nwf) { // wavefront loop (x direction)
            if(xe-xi <= nwf){
                nwf = xe-xi;
                end =1;
            }
            yb = yb_r;
            ye = ye_r;
            kt = xi;
            int kte=kt+nwf;

            float * __restrict v1_v;
            float * __restrict v2_v;
            float * __restrict v3_v;
            float * __restrict u1_v;
            float * __restrict u2_v;
            float * __restrict u3_v;
            const float * __restrict coef0_v;
            int t_real=0;
            int tb_real=(tb)/2+1;
            int t0_real=(t0)/2+1;
            for(int t=tb; t< te; t++){ // Diamond blocking in time
                t_real=(t)/2+1;
                hFloat* output_buffer = NULL;  //@KADIR
                int mod = (t)%2;
//                MSG("t=%d",t);
                if(mod){ // compute v from p
                    u1 = p11 ; //p
                    u2 = p12 ;
                    u3 = p13 ;
                    v1 = p21 ; //vx
                    v2 = p22 ; //vy
                    v3 = p23 ; //vz
//#pragma omp barrier
                    const Myfloat coef=stencil_ctx.dt;
                    for(int ix=kt; ix<kte; ix++){    // X
                        if( ((ix)/th_nwf)%th_x == tid_x ) {
                            for(int iy=yb; iy<ye; iy++) {
                                {
                                    v1_v = &(v1[ix*nnyz+iy*nnz]);
                                    v2_v = &(v2[ix*nnyz+iy*nnz]);
                                    v3_v = &(v3[ix*nnyz+iy*nnz]);
                                    u1_v = &(u1[ix*nnyz+iy*nnz]);
                                    u2_v = &(u2[ix*nnyz+iy*nnz]);
                                    u3_v = &(u3[ix*nnyz+iy*nnz]);
#pragma ivdep
                                    for(int iz=ib; iz<ie; iz++) {
                                        const Myfloat xum4 = u1_v[-3*nnyz + iz];
                                        const Myfloat xum3 = u1_v[-2*nnyz + iz];
                                        const Myfloat xum2 = u1_v[-1*nnyz + iz];
                                        const Myfloat xum1 = u1_v[ 0*nnyz + iz];
                                        const Myfloat xu0  = u1_v[ 1*nnyz + iz];
                                        const Myfloat xup1 = u1_v[ 2*nnyz + iz];
                                        const Myfloat xup2 = u1_v[ 3*nnyz + iz];
                                        const Myfloat xup3 = u1_v[ 4*nnyz + iz];

                                        Myfloat d_pr_x  = ( ( FDM_O1_8_2_A1 * (xu0  - xum1)
                                                              + FDM_O1_8_2_A2 * (xup1 - xum2)
                                                              + FDM_O1_8_2_A3 * (xup2 - xum3)
                                                              + FDM_O1_8_2_A4 * (xup3 - xum4)) * inv_dx) ;

                                        v1_v[iz] += stencil_ctx.dt * d_pr_x;

                                        const Myfloat yum4 = u1_v[-3*nnz + iz];
                                        const Myfloat yum3 = u1_v[-2*nnz + iz];
                                        const Myfloat yum2 = u1_v[-1*nnz + iz];
                                        const Myfloat yum1 = u1_v[ 0*nnz + iz];
                                        const Myfloat yu0  = u1_v[ 1*nnz + iz];
                                        const Myfloat yup1 = u1_v[ 2*nnz + iz];
                                        const Myfloat yup2 = u1_v[ 3*nnz + iz];
                                        const Myfloat yup3 = u1_v[ 4*nnz + iz];

                                        Myfloat d_pr_y  = ( ( FDM_O1_8_2_A1 * (yu0  - yum1)
                                                              + FDM_O1_8_2_A2 * (yup1 - yum2)
                                                              + FDM_O1_8_2_A3 * (yup2 - yum3)
                                                              + FDM_O1_8_2_A4 * (yup3 - yum4)) * inv_dy) ;

                                        v2_v[iz] += stencil_ctx.dt * d_pr_y;

                                        const Myfloat zum4 = u1_v[-3 + iz];
                                        const Myfloat zum3 = u1_v[-2 + iz];
                                        const Myfloat zum2 = u1_v[-1 + iz];
                                        const Myfloat zum1 = u1_v[ 0 + iz];
                                        const Myfloat zu0  = u1_v[ 1 + iz];
                                        const Myfloat zup1 = u1_v[ 2 + iz];
                                        const Myfloat zup2 = u1_v[ 3 + iz];
                                        const Myfloat zup3 = u1_v[ 4 + iz];

                                        Myfloat d_pr_z  = ( ( FDM_O1_8_2_A1 * (zu0  - zum1)
                                                              + FDM_O1_8_2_A2 * (zup1 - zum2)
                                                              + FDM_O1_8_2_A3 * (zup2 - zum3)
                                                              + FDM_O1_8_2_A4 * (zup3 - zum4)) * inv_dz) ;

                                        v3_v[iz] += stencil_ctx.dt * d_pr_z;
                                    }
                                }
                            }

                        }
                    }
                } else{// compute p from v
                    u1=	p21 ; //vx
                    u2=	p22 ; //vy
                    u3=	p23 ; //vz
                    v1=	p11 ; //p
                    v2=	p12 ;
                    v3=	p13 ;
                    for(int ix=kt; ix<kte; ix++){
                        if( ((ix)/th_nwf)%th_x == tid_x ) {
                            for(int iy=yb; iy<ye; iy++) {
//                                printf("iy=%d\n",iy);
                                v1_v = &(v1[ix*nnyz+iy*nnz]);
                                v2_v = &(v2[ix*nnyz+iy*nnz]);
                                v3_v = &(v3[ix*nnyz+iy*nnz]);
                                u1_v = &(u1[ix*nnyz+iy*nnz]);
                                u2_v = &(u2[ix*nnyz+iy*nnz]);
                                u3_v = &(u3[ix*nnyz+iy*nnz]);
//                                coef0_v = &(roc2[ix*nnyz+iy*nnz]); original
                                coef0_v = &(roc2[(ix-NHALO)*nnyz_v+(iy-NHALO)*nnz_v]);
#pragma ivdep
                                for(int iz=ib; iz<ie; iz++) {
                                    const Myfloat xum4 = u1_v[-4*nnyz + iz];
                                    const Myfloat xum3 = u1_v[-3*nnyz + iz];
                                    const Myfloat xum2 = u1_v[-2*nnyz + iz];
                                    const Myfloat xum1 = u1_v[-1*nnyz + iz];
                                    const Myfloat xu0  = u1_v[ 0*nnyz + iz];
                                    const Myfloat xup1 = u1_v[ 1*nnyz + iz];
                                    const Myfloat xup2 = u1_v[ 2*nnyz + iz];
                                    const Myfloat xup3 = u1_v[ 3*nnyz + iz];

                                    Myfloat d_vx_x  = ( ( FDM_O1_8_2_A1 * (xu0  - xum1)
                                                          + FDM_O1_8_2_A2 * (xup1 - xum2)
                                                          + FDM_O1_8_2_A3 * (xup2 - xum3)
                                                          + FDM_O1_8_2_A4 * (xup3 - xum4)) * inv_dx) ;

                                    const Myfloat yum4 = u2_v[-4*nnz + iz];
                                    const Myfloat yum3 = u2_v[-3*nnz + iz];
                                    const Myfloat yum2 = u2_v[-2*nnz + iz];
                                    const Myfloat yum1 = u2_v[-1*nnz + iz];
                                    const Myfloat yu0  = u2_v[ 0*nnz + iz];
                                    const Myfloat yup1 = u2_v[ 1*nnz + iz];
                                    const Myfloat yup2 = u2_v[ 2*nnz + iz];
                                    const Myfloat yup3 = u2_v[ 3*nnz + iz];

                                    Myfloat d_vy_y  = ( ( FDM_O1_8_2_A1 * (yu0  - yum1)
                                                          + FDM_O1_8_2_A2 * (yup1 - yum2)
                                                          + FDM_O1_8_2_A3 * (yup2 - yum3)
                                                          + FDM_O1_8_2_A4 * (yup3 - yum4)) * inv_dy) ;

                                    const Myfloat zum4 = u3_v[-4 + iz];
                                    const Myfloat zum3 = u3_v[-3 + iz];
                                    const Myfloat zum2 = u3_v[-2 + iz];
                                    const Myfloat zum1 = u3_v[-1 + iz];
                                    const Myfloat zu0  = u3_v[ 0 + iz];
                                    const Myfloat zup1 = u3_v[ 1 + iz];
                                    const Myfloat zup2 = u3_v[ 2 + iz];
                                    const Myfloat zup3 = u3_v[ 3 + iz];

                                    Myfloat d_vz_z  = ( ( FDM_O1_8_2_A1 * (zu0  - zum1)
                                                          + FDM_O1_8_2_A2 * (zup1 - zum2)
                                                          + FDM_O1_8_2_A3 * (zup2 - zum3)
                                                          + FDM_O1_8_2_A4 * (zup3 - zum4)) * inv_dz);
                                    v3_v[iz] += coef0_v[iz] * (d_vx_x + d_vy_y + d_vz_z);
                                    v3_v[iz] = dampx[ix+lstencil] * v3_v[iz];
                                    v3_v[iz] = dampy[iy+lstencil] * v3_v[iz];
                                    v3_v[iz] = dampz[iz+lstencil] * v3_v[iz];
                                }
                                ///////  save sismos
                                //////////////////////////////////////////
//                                if (data->rcv_len>0){
//                                    data->sismos[data->rcv_len*(t0_real+(t_real-tb_real))+(iy-4)*(nnx-2*lstencil)+(ix-4)]=(v3_v[iz_]);
//                                }
                                //////////////////////////////////////////
                            }
                            ///
                            if( (gp->source_point_enabled==1)
                                && (gp->lsource_pt[2] >= ib ) //@KADIR
                                && (gp->lsource_pt[2] <  ie ) //@KADIR
                                && (gp->lsource_pt[1] >= yb ) //@KADIR
                                && (gp->lsource_pt[1] <  ye ) //@KADIR
                                && (gp->lsource_pt[0] == ix ) //@KADIR
                                    )
                            {
//								gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))] = F2H(H2F(gp->U1[((1ULL)*((gp->lsource_pt[2])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[0]))]) +gp->src_exc_coef[isrc_exc]);//@KADIR
                                gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))] = F2H(H2F(gp->U1[((1ULL)*((gp->lsource_pt[0])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[2]))]) + gp->src_exc_coef[isrc_exc]);//@KADIR

                                if(0)  printf("DIA\tts:%d idzU:-- valU:%.4f src_exc_coef:%.4f coef:%g %g %g %g %g\ti(%d-%d) %d/%d\n", isrc_exc, H2F(gp->U1[((1ULL)*((gp->lsource_pt[2])*(gp->ldomain_shape[1])+( gp->lsource_pt[1]))*(gp->ldomain_shape[0])+(gp->lsource_pt[0]))]),  gp->src_exc_coef[isrc_exc], coef[0], coef[1], coef[2], coef[3], coef[4],
                                              ib, ie, omp_get_thread_num(), omp_get_num_threads());
                                isrc_exc++;
                            }
                        }
                    }
                }


                // Update block size in Y
                if(t< t_dim){ // lower half of the diamond
                    yb += -b_inc;
                    ye += e_inc;
                }else{ // upper half of the diamond
                    yb += b_inc;
                    ye += -e_inc;
                }
                kte=max(kte-NHALO,xb);
                if (end==1) kte =xe;
                kt=max(kt-NHALO,xb);
                t_start = get_wall_time();
#pragma omp barrier
                stencil_ctx.t_wait[gtid] += get_wall_time() - t_start;
            } // diamond blocking in time (time loop)
        } // wavefront loop
    } // parallel region
}
*/
