Here is an example of how to validate that the GPU is generating the same output as the CPU:
```
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave$ ./bin/modeling 
[SIMWAVE MSG]:... dealing with shot number 0 (source @ 4495).
[SIMWAVE MSG]:... dealing with shot number 1 (source @ 4545).
[SIMWAVE MSG]:... dealing with shot number 2 (source @ 5935).
[SIMWAVE MSG]:... dealing with shot number 3 (source @ 5985).
[SIMWAVE MSG]:... dealing with shot number 4 (source @ 7375).
[SIMWAVE MSG]:... dealing with shot number 5 (source @ 7425).
[SIMWAVE MSG]:... dealing with shot number 6 (source @ 8815).
[SIMWAVE MSG]:... dealing with shot number 7 (source @ 8865).
[SIMWAVE MSG]:... dealing with shot number 8 (source @ 10255).
[SIMWAVE MSG]:... dealing with shot number 9 (source @ 10305).
[SIMWAVE MSG]:... dealing with shot number 10 (source @ 11695).
[SIMWAVE MSG]:... dealing with shot number 11 (source @ 11745).
[SIMWAVE MSG]:... dealing with shot number 12 (source @ 13135).
[SIMWAVE MSG]:... dealing with shot number 13 (source @ 13185).
[SIMWAVE MSG]:... dealing with shot number 14 (source @ 14575).
[SIMWAVE MSG]:... dealing with shot number 15 (source @ 14625).
[SIMWAVE MSG]:... dealing with shot number 16 (source @ 16015).
[SIMWAVE MSG]:... dealing with shot number 17 (source @ 16065).
[SIMWAVE MSG]:... dealing with shot number 18 (source @ 17455).
[SIMWAVE MSG]:... dealing with shot number 19 (source @ 17505).
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave$ ./bin/modeling --cpu 
[SIMWAVE MSG]:... dealing with shot number 0 (source @ 4495).
[SIMWAVE MSG]:... dealing with shot number 1 (source @ 4545).
[SIMWAVE MSG]:... dealing with shot number 2 (source @ 5935).
[SIMWAVE MSG]:... dealing with shot number 3 (source @ 5985).
[SIMWAVE MSG]:... dealing with shot number 4 (source @ 7375).
[SIMWAVE MSG]:... dealing with shot number 5 (source @ 7425).
[SIMWAVE MSG]:... dealing with shot number 6 (source @ 8815).
[SIMWAVE MSG]:... dealing with shot number 7 (source @ 8865).
[SIMWAVE MSG]:... dealing with shot number 8 (source @ 10255).
[SIMWAVE MSG]:... dealing with shot number 9 (source @ 10305).
[SIMWAVE MSG]:... dealing with shot number 10 (source @ 11695).
[SIMWAVE MSG]:... dealing with shot number 11 (source @ 11745).
[SIMWAVE MSG]:... dealing with shot number 12 (source @ 13135).
[SIMWAVE MSG]:... dealing with shot number 13 (source @ 13185).
[SIMWAVE MSG]:... dealing with shot number 14 (source @ 14575).
[SIMWAVE MSG]:... dealing with shot number 15 (source @ 14625).
[SIMWAVE MSG]:... dealing with shot number 16 (source @ 16015).
[SIMWAVE MSG]:... dealing with shot number 17 (source @ 16065).
[SIMWAVE MSG]:... dealing with shot number 18 (source @ 17455).
[SIMWAVE MSG]:... dealing with shot number 19 (source @ 17505).
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave$ ls data/
augmented_vel.raw  gpu_sismos_13.raw  gpu_sismos_18.raw  gpu_sismos_4.raw  gpu_sismos_9.raw  sismos_12.raw  sismos_17.raw  sismos_3.raw  sismos_8.raw
gpu_sismos_0.raw   gpu_sismos_14.raw  gpu_sismos_19.raw  gpu_sismos_5.raw  pml.raw           sismos_13.raw  sismos_18.raw  sismos_4.raw  sismos_9.raw
gpu_sismos_10.raw  gpu_sismos_15.raw  gpu_sismos_1.raw   gpu_sismos_6.raw  sismos_0.raw      sismos_14.raw  sismos_19.raw  sismos_5.raw  source.raw
gpu_sismos_11.raw  gpu_sismos_16.raw  gpu_sismos_2.raw   gpu_sismos_7.raw  sismos_10.raw     sismos_15.raw  sismos_1.raw   sismos_6.raw  source.txt
gpu_sismos_12.raw  gpu_sismos_17.raw  gpu_sismos_3.raw   gpu_sismos_8.raw  sismos_11.raw     sismos_16.raw  sismos_2.raw   sismos_7.raw
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave$ cd scripts/
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave/scripts$ make
gcc -g diff_to.c -o diff_to -lm
gcc -g isnan.c -o isnan -lm
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave/scripts$ cd ..
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave$ ./scripts/isnan data/sismos_0.raw 
...
... file: data/sismos_0.raw
... size: 400000
... N   : 100000
... min : -1.23975 @ 10700
... max : 1.3109 @ 8800
... rms : 0.0238944
... nan : 0 (0.000000 %)
...
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave$ ./scripts/isnan data/gpu_sismos_0.raw 
...
... file: data/gpu_sismos_0.raw
... size: 400000
... N   : 100000
... min : -1.23974 @ 10700
... max : 1.3109 @ 8800
... rms : 0.0238945
... nan : 0 (0.000000 %)
...
isaid@nvparis-isaid-dt:~/code/gitlab/nvidia-kaust/simwave$ ./scripts/diff_to data/gpu_sismos_0.raw data/sismos_0.raw 
...
... [DIFF TAB] generated is : data/diff_gpu_sismos_0.raw (100000 elements)
...
... [GPU_TAB] min:            -1.239745 @ [               10700] | max:             1.310899 @ [                8800] | rms:           0.02389446
... [REF_TAB] min:            -1.239746 @ [               10700] | max:             1.310897 @ [                8800] | rms:           0.02389445
... [DIF_TAB] min:        -1.549721e-06 @ [                8900] | max:         2.026558e-06 @ [               11300] | rms:          4.20505e-08
...
```
