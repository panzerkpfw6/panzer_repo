##
## @copyright Copyright 2017- Issam Said. All rights reserved.
## This file is part of \b stencil.
##
## \b stencil is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## stencil is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with \b stencil.  If not, see <http://www.gnu.org/licenses/>.
##
## @file CMakeLists.txt
## @author Issam SAID
## @brief Install all the subdirectories when make install is invoked.
##

# ## Install the hiCL Fortran interface.
# if (EXISTS ${HICL_FORTRAN_BUILD_LIB})
#     message(STATUS "Post-install: install the Fortran interface")
#     file(COPY ${HICL_FORTRAN_BUILD_LIB}
#          DESTINATION ${HICL_FORTRAN_INSTALL_DIR})
# endif (EXISTS ${HICL_FORTRAN_BUILD_LIB})

# ## Install the hiCL C unit testing executable.
# if (EXISTS ${HICL_TEST_BUILD_EXE})
#     message(STATUS "Post-install: install the hiCL C unit test binary")
#     file(COPY ${HICL_TEST_BUILD_EXE}
#          DESTINATION ${HICL_TEST_INSTALL_DIR})
# endif (EXISTS ${HICL_TEST_BUILD_EXE})

# ## Install the hiCL Fortran interface unit testing executable.
# if (EXISTS ${HICL_TEST_FORTRAN_BUILD_EXE})
#     message(STATUS "Post-install: install the hiCL Fortran unit test binary")
#     file(COPY ${HICL_TEST_FORTRAN_BUILD_EXE}
#          DESTINATION ${HICL_TEST_FORTRAN_INSTALL_DIR})
# endif (EXISTS ${HICL_TEST_FORTRAN_BUILD_EXE})

# ## Install the hiCL C examples.
# if (DEFINED HICL_C_EXAMPLES_FILES)
#     foreach (file ${HICL_C_EXAMPLES_FILES})
#         if (EXISTS ${file})
#             get_filename_component(n ${file} NAME)
#             message(STATUS "Post-install: install the hiCL C example ${n}")
#             file(COPY ${file} 
#                  DESTINATION ${HICL_C_EXAMPLES_INSTALL_DIR})
#         endif (EXISTS ${file})
#     endforeach(file ${HICL_C_EXAMPLES_FILES})
# endif (DEFINED HICL_C_EXAMPLES_FILES)

# ## Install the hiCL Fortran examples.
# if (DEFINED HICL_FORTRAN_EXAMPLES_FILES)
#     foreach (file ${HICL_FORTRAN_EXAMPLES_FILES})
#         if (EXISTS ${file})
#             get_filename_component(n ${file} NAME)
#             message(STATUS "Post-install: install the hiCL Fortran example ${n}")
#             file(COPY ${file} 
#                  DESTINATION ${HICL_FORTRAN_EXAMPLES_INSTALL_DIR})
#         endif (EXISTS ${file})
#     endforeach(file ${HICL_FORTRAN_EXAMPLES_FILES})
# endif (DEFINED HICL_FORTRAN_EXAMPLES_FILES)