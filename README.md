```
This branch represents acoustic code based on inhomogeneous density with XYZ ordering imported from stencil-dev code.

NAME
	stencil - a synthetic wave simulator
SYNTAX
	stencil [options]

OPTIONS
	-f, --file         VALUE       parse parameters from file.
	-v, --verbose                  enable verbose mode.
	-c, --cpu                      run serial code on the CPU.
	-l, --local        V0,V1,...   set the GPU block dimensions.
	-d, --device       VALUE       select the GPU device.
	-o, --one                      use only one GPU kernel.
	    --gpu_options  VALUE       set GPU build options.
	-m, --ms                       use milliseconds for timing.
	-u, --us                       use microseconds for timing.
	-n, --ns                       use nanoseconds  for timing.
	    --grid         V0,V1,...   domain x dimension.
	    --dgrid        V0,V1,...   space delta step.
	    --source_loc   V0,V1,...   source location in the grid.
	-i, --iter         VALUE       simulation time step number.
	    --cfl          VALUE       CFL percentage.
	    --fmax         VALUE       source max frequency.
	    --vmin         VALUE       min velocity.
	    --vmax         VALUE       max velocity.
	    --output       VALUE       snapshots file path.
	    --nosnap                   do not save snapshots.
	    --nbsnap       VALUE       snapshot frequency.
	    --check                    check the GPU results.
	    --taper        V0,V1,...   PML tapers.
	-e, --epsilon      VALUE       the margin of floating point errors.
We are compiling the project with cmake software.
#### CMAKE requirements
Please edit CMAKE_SOURCE_DIR,CMAKE_BINARY_DIR filepaths in ./Makefile  
Also edit CMakeLists.txt files in ./, ./app and ./src subdirectories:
./CMakeLists.txt
./app/CMakeLists.txt
./src/CMakeLists.txt
```
