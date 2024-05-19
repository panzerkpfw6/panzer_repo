## supported CUDA architectures.
set (SM60 -gencode=arch=compute_60,code=sm_60)
set (SM61 -gencode=arch=compute_61,code=sm_61)
set (SM62 -gencode=arch=compute_62,code=sm_62)
set (SM70 -gencode=arch=compute_70,code=sm_70)
set (SM75 -gencode=arch=compute_75,code=sm_75)

## if compiling against CUDA Toolkit 8.x
if(CUDA_VERSION_MAJOR MATCHES 8)
  if(CUDA_VERSION_MINOR MATCHES 0)
    set(CUDA_ARCH "SM60 SM61 SM62" 
        CACHE STRING
        "target arch: SM60 SM61 SM62")
  else(CUDA_VERSION_MINOR MATCHES 0)  
    set(CUDA_ARCH "SM60 SM61 SM62" 
        CACHE STRING
        "target arch: SM60 SM61 SM62")
  endif(CUDA_VERSION_MINOR MATCHES 0)  
endif(CUDA_VERSION_MAJOR MATCHES 8)

## if compiling against CUDA Toolkit 9.x
if(CUDA_VERSION_MAJOR MATCHES 9)
  if(CUDA_VERSION_MINOR MATCHES 0)
    set(CUDA_ARCH "SM60 SM61 SM62 SM70" 
        CACHE STRING
        "target arch: SM60 SM61 SM62 SM70")
  else(CUDA_VERSION_MINOR MATCHES 0)  
    set(CUDA_ARCH "SM60 SM61 SM62 SM70" 
        CACHE STRING
        "target arch: SM60 SM61 SM62 SM70")
  endif(CUDA_VERSION_MINOR MATCHES 0)  
endif(CUDA_VERSION_MAJOR MATCHES 9)

## if compiling against CUDA Toolkit 10.x
if(CUDA_VERSION_MAJOR MATCHES 10)
  if(CUDA_VERSION_MINOR MATCHES 0)
    set(CUDA_ARCH "SM60 SM70 SM75" 
        CACHE STRING
        "target arch: SM60 SM70 SM75")
  else(CUDA_VERSION_MINOR MATCHES 0)  
    set(CUDA_ARCH "SM60 SM70 SM75" 
        CACHE STRING
        "target arch: SM60 SM70 SM75")
  endif(CUDA_VERSION_MINOR MATCHES 0)  
endif(CUDA_VERSION_MAJOR MATCHES 10)

## replace ' ' with ; to match the proper cmake format
string(REGEX REPLACE " " ";" CUDA_ARCH ${CUDA_ARCH})

## set the compiler flags for each NV target
foreach(target ${CUDA_ARCH})
	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${${target}}")
endforeach(target ${CUDA_ARCH})

## set the host compiler if any
if (CMAKE_CUDA_HOST_COMPILER)
	set(CMAKE_CUDA_FLAGS "-ccbin=${CMAKE_CUDA_HOST_COMPILER} ${CMAKE_CUDA_FLAGS}")
endif (CMAKE_CUDA_HOST_COMPILER)

