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
                    ./SB_1st.out 512 512 512 10 100 0 0
                    ## simulate one shot in the center of domain
                    nx=128;ny=256;nz=512;
                    nt=10;  dt=0.001;
                    export shot=16447;  # position of the source in x,y coordinates.check ./data/acquisition.txt
                    export src_depth=256;
                    ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2 --dshot 1 --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 1 --fmax 8;
                    ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  --dshot 1 --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 2 --fmax 8;
                    '''
            }
        }
        stage ('test_TB') {
            steps {
                sh '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5;
                    nx=128;ny=256;nz=512;
                    nt=57; dt=0.001;
                    x=2 y=2 z=1 t=7 w=20;
                    tgs=$(expr $x * $y * $z);echo $tgs;
                    ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $nt --tgs $tgs --nb_thg $(expr $OMP_NUM_THREADS / $tgs) --thx $x --thy $y --thz $z --tdim $t --nwf $w --mode 2  --dshot 1 --first $shot --last $shot  --fwd_steps 3 -c --order 1 --fmax 8;
                    ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $nt --tgs $tgs --nb_thg $(expr $OMP_NUM_THREADS / $tgs) --thx $x --thy $y --thz $z --tdim $t --nwf $w --mode 2  --dshot 1 --first $shot --last $shot  --fwd_steps 3 -c --order 2 --fmax 8;
                    '''
            }
        }
        stage ('compare_wavefields_for_SB_TB') {
            steps {
                sh '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    export OMP_NUM_THREADS=4
                    ####################################################
                    echo "Compare wavefields"
                    export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    export TIME_SB_1st=504 #@pavel in SB the nt should one time less than in correponding TB.

                    ### grid size 512*512*512
                    x=2 y=2 z=1 t=7 w=20;  numactl --interleave=all ./TB_1st.out 512 512 512 $TIME_TB_1st 1 $x $y $z $t $w 0 1;    # we record last time sample wavefield now and set is to 1
                    ./SB_1st.out 512 512 512 $TIME_SB_1st $TIME_SB_1st 0 1;
                    ./scripts/diff_to ./TB_1st.raw ./SB_1st.raw
                    ### grid size 256*256*256
                    x=2 y=2 z=1 t=7 w=20;  numactl --interleave=all ./TB_1st.out 256 256 256 $TIME_TB_1st 1 $x $y $z $t $w 0 1;    # we record last time sample wavefield now and set is to 1
                    ./SB_1st.out 256 256 256 $TIME_SB_1st $TIME_SB_1st 0 1;
                    ./scripts/diff_to ./TB_1st.raw ./SB_1st.raw
                    '''
            }
        }
        stage ('compare_wavefields_for_SBabc_TBabc') {
            steps {
                sh '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    export OMP_NUM_THREADS=4
                    ####################################################
                    export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    export TIME_SB_1st=504 #@pavel in SB the nt should one time less than in correponding TB.
                    ### grid size 256*256*256
                    x=2 y=2 z=1 t=7 w=20;
                    ./SB_1st-abc.out 256 256 256 $TIME_SB_1st $TIME_SB_1st 0 1 0 0;
                    numactl --interleave=all ./TB_1st-abc.out 256 256 256 $TIME_TB_1st 1 $x $y $z $t $w 0 1 0 0;
                    ./scripts/diff_to ./TB_1st_abc.raw ./SB_1st_abc.raw
                    '''
            }
        }
        stage ('compare_wavefields_for_SB_TB_2nd_order') {
            steps {
                sh  '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    export OMP_NUM_THREADS=4
                    ####################################################
                    export TIME_TB_2nd=514 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    export TIME_SB_2nd=513 #@pavel in SB the nt should one time less than in correponding TB.
                    ### grid size 256*256*256
                    x=2 y=2 z=1 t=7 w=20;  numactl --interleave=all ./TB_2nd.out 256 256 256 $TIME_TB_2nd 1 $x $y $z $t $w 0 1;    # we record last time sample wavefield now and set is to 1
                    ./SB_2nd.out 256 256 256 $TIME_SB_2nd $TIME_SB_2nd 0 1;
                    echo 'compare_wavefields_for_SB_TB_2nd_order'
                    ./scripts/diff_to ./TB_2nd.raw ./SB_2nd.raw
                    '''
            }
        }
        stage ('compare_wavefields_for_SB_TB_2nd_order_abc') {
            steps {
                sh  '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    export OMP_NUM_THREADS=4
                    ####################################################
                    export TIME_TB_2nd=514 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    export TIME_SB_2nd=513 #@pavel in SB the nt should one time less than in correponding TB.
                    ### grid size 256*256*256
                    x=2 y=2 z=1 t=7 w=20;  numactl --interleave=all ./TB_2nd-abc.out 256 256 256 $TIME_TB_2nd 1 $x $y $z $t $w 0 1;    # we record last time sample wavefield now and set is to 1
                    ./SB_2nd-abc.out 256 256 256 $TIME_SB_2nd $TIME_SB_2nd 0 1;
                    echo 'compare_wavefields_for_SB_TB_2nd_order_abc'
                    ./scripts/diff_to ./TB_2nd_abc.raw ./SB_2nd_abc.raw
                    '''
            }
        }
        stage ('test_sismos_options_for_SB') {
            steps {
                sh  '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    export OMP_NUM_THREADS=4
                    ####################################################
                    export TIME_SB_1st=54 #@pavel in SB the nt should one time less than in correponding TB.
                    ### we test functionality of 2 bool flags "read_v,shot_flag". Does the program crash or not.
                    # first two flags are write_velocity_grid,write_snapshot
                    ./SB_1st-abc.out 256 256 256 $TIME_SB_1st $TIME_SB_1st 0 1 1 0;
                    ./SB_1st-abc.out 256 256 256 $TIME_SB_1st $TIME_SB_1st 0 1 1 1;
                    ./SB_1st-abc.out 256 256 256 $TIME_SB_1st $TIME_SB_1st 0 1 0 0;
                    ./SB_1st-abc.out 256 256 256 $TIME_SB_1st $TIME_SB_1st 0 1 0 1;
                    '''
            }
        }
        stage ('test_sismos_options_for_TB') {
            steps {
                sh  '''#!/bin/bash -le
                    ####################################################
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    export OMP_NUM_THREADS=4
                    ####################################################
                    export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    ### we test functionality of 2 bool flags "read_v,shot_flag". Does the program crash or not.
                    # first two flags are write_velocity_grid,write_snapshot
                    x=2 y=2 z=1 t=7 w=20;
                    echo '!!!!!!!0 0!!!!!!!'
                    numactl --interleave=all  ./TB_1st-abc.out 256 256 256 $TIME_TB_1st 1 $x $y $z $t $w 0 1 0 0;
                    echo '!!!!!!!0 1!!!!!!!'
                    numactl --interleave=all ./TB_1st-abc.out 256 256 256 $TIME_TB_1st 1 $x $y $z $t  $w 0 1 0 1; #ok 505
                    echo '!!!!!!!1 0!!!!!!!'
                    numactl --interleave=all ./TB_1st-abc.out 256 256 256 $TIME_TB_1st 1 $x $y $z $t $w 0 1 1 0;
                    echo '!!!!!!! 1 1 !!!!!!!'
                    numactl --interleave=all ./TB_1st-abc.out 256 256 256 $TIME_TB_1st 1 $x $y $z $t $w 0 1 1 1; #ok 505
                    '''
            }
        }
        stage ('compare_sismos_for_SBabc_TBabc') {
            steps {
                sh  '''#!/bin/bash -le
                    module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
                    export OMP_NUM_THREADS=4
                    ####################################################
                    export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
                    export TIME_SB_1st=504 #@pavel in SB the nt should one time less than in correponding TB.
                    ####################################################
                    #######   test with read_v=0   #######
                    cd ./SB_abc;
                    pwd
                    icpc -xHost -qopenmp -O3 -I. test_SB_kernel.cpp -o ./SB_1st-abc.out
                    ./SB_1st-abc.out 256 512 256 $TIME_SB_1st $TIME_SB_1st 0 0 0 1;

                    cd ../TB_abc;
                    icpc -xHost -qopenmp -O3 -I. test_TB_kernel.cpp -o ./TB_1st-abc.out
                    x=2 y=2 z=1 t=7 w=20;
                    ./TB_1st-abc.out 256 512 256 $TIME_TB_1st 1 $x $y $z $t $w 0 0 0 1;
                    cd ..

                    python test_tb_sismos.py './SB_abc/data/sismos.raw' './TB_abc/data/sismos.raw' 256 256 504 505
                    #######   test with read_v=1   #######
                    cd ./SB_abc;
                    icpc -xHost -qopenmp -O3 -I. test_SB_kernel.cpp -o ./SB_1st-abc.out
                    ./SB_1st-abc.out 256 512 256 $TIME_SB_1st $TIME_SB_1st 0 0 1 1;

                    cd ../TB_abc;
                    icpc -xHost -qopenmp -O3 -I. test_TB_kernel.cpp -o ./TB_1st-abc.out
                    x=2 y=2 z=1 t=7 w=20;
                    ./TB_1st-abc.out 256 512 256 $TIME_TB_1st 1 $x $y $z $t $w 0 0 1 1;
                    cd ..

                    python test_tb_sismos.py './SB_abc/data/sismos.raw' './TB_abc/data/sismos.raw' 256 512 504 505
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
