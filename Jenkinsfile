pipeline {
    agent { label 'jenkinsfile' }
    triggers {
        pollSCM('H/10 * * * *')
    }

    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }

    stages {
        stage ('compile_project') {
            steps {
                sh '''#!/bin/bash -le
                    exec > >(while read line; do echo "$(date): $line"; done | tee log-modeling_.log) 2>&1
                    echo $hostname
                    export OMP_NUM_THREADS=4
                    export OMP_PROC_BIND=true
                    export OMP_PLACES=threads
                    export OMP_NESTED='True'
                    export granularity=fine
                    export KMP_AFFINITY=compact
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    module load cmake

                    pwd
                    ls
                    mv -f ./CMakeCache.txt ./CMakeCache-old.txt
                    CC=icc CXX=icpc cmake .
                    make clean
                    make VERBOSE=1
                    make install
                '''
            }
        }
        stage ('test_SB') {
            steps {
                sh '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    ## simulate one shot in the center of domain
                    nx=128;ny=256;nz=512;
                    nt=10;  dt=0.001;
                    export shot=16447;  # position of the source in x,y coordinates.check ./data/acquisition.txt
                    export src_depth=256;
                    ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --mode 2 --dshot 1 --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 1 --fmax 8;
                    ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --mode 2  --dshot 1 --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 2 --fmax 8;
                    '''
            }
        }
        stage ('test_TB') {
            steps {
                sh '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5;
                    nx=128;ny=256;nz=512;
                    nt=57; dt=0.001;
                    x=2; y=2; z=1; t=7; w=20; tgs=4;
                    export OMP_NUM_THREADS=4
                    export shot=16447;  # position of the source in x,y coordinates.check ./data/acquisition.txt
                    export src_depth=256;
                    ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $nt --tb_thread_group_size $tgs \
                     --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
                     --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c --src_depth $src_depth --order 1 --fmax 8;
                    ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $nt --tb_thread_group_size $tgs \
                     --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
                     --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c --src_depth $src_depth --order 2 --fmax 8;
                    '''
            }
        }
        stage ('compare_wavefields_for_SB_TB') {
            steps {
                sh '''#!/bin/bash -le
                    echo "compare_wavefields_for_SB_TB"
                    export OMP_NUM_THREADS=4 #4
                    export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    export TIME_SB_1st=505 #@pavel in SB the nt should one time less than in correponding TB.
                    ### grid size 256*256*256
                    nx=256;ny=256;nz=256;
                    x=2; y=2; z=1; t=7; w=20; tgs=4;
                    export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
                    export src_depth=128;
                    ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
                     --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
                     --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
                     --src_depth $src_depth --order 1 --fmax 8;
                    ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2 --dshot 1 \
                      --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 1 --fmax 8;
                    ./scripts_useful/diff_to ./snapshot_TB1st_505 ./snapshot_SB1st_505
                    '''
            }
        }
        stage ('compare_wavefields_for_SB_TB_2nd_order') {
            steps {
                sh  '''#!/bin/bash -le
                    echo "compare_wavefields_for_SB_TB_2nd_order"
                    export OMP_NUM_THREADS=4
                    export TIME_TB_2nd=514 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    export TIME_SB_2nd=514 #@pavel in SB the nt should one time less than in correponding TB.
                    ### grid size 256*256*256
                    nx=256;ny=256;nz=256;
                    x=2; y=2; z=1; t=7; w=20; tgs=4;
                    export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
                    export src_depth=128;
                    ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_2nd --tb_thread_group_size $tgs \
                     --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
                     --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
                     --src_depth $src_depth --order 2 --fmax 8;
                    ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_2nd --mode 2 --dshot 1 \
                      --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 2 --fmax 8;
                    ./scripts_useful/diff_to ./snapshot_TB2nd_514 ./snapshot_SB2nd_514
                    '''
            }
        }
        stage ('compare_wavefields_for_SB_1st_2nd_order') {
            steps {
                sh  '''#!/bin/bash -le
                    echo "compare_wavefields_for_SB_1st_2nd_order"
                    export OMP_NUM_THREADS=4
                    export TIME_SB_1st=520 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    export TIME_SB_2nd=520 #@pavel in SB the nt should one time less than in correponding TB.
                    ### grid size 256*256*256
                    nx=256;ny=256;nz=256;
                    export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
                    export src_depth=128;
                    ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2 --dshot 1 \
                      --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 1 --fmax 8;
                    ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_2nd --mode 2 --dshot 1 \
                      --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 2 --fmax 8;
                    ./scripts_useful/diff_to ./snapshot_SB1st_520 ./snapshot_SB2nd_520
                    '''
            }
        }
        stage ('test_sismos_options_for_SB') {
            steps {
                sh  '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    '''
            }
        }
        stage ('test_sismos_options_for_TB') {
            steps {
                sh  '''#!/bin/bash -le
                    '''
            }
        }
        stage ('compare_sismos_for_SBabc_TBabc') {
            steps {
                sh  '''#!/bin/bash -le
                    '''
            }
        }
    }
    // Post build actions
    post {
        //always {
        //}
        //success {
        //}
        //unstable {
        //}
        //failure {
        //}
        unstable {
                emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build is UNSTABLE", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
        }
        failure {
                emailext body: "${env.JOB_NAME} - Please go to ${env.BUILD_URL}", subject: "Jenkins Pipeline build FAILED", recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'RequesterRecipientProvider']]
        }
    }
}
