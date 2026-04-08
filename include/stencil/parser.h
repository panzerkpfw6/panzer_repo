#ifndef __STENCIL_PARSER_H_
#define __STENCIL_PARSER_H_

///
/// @copyright Copyright 2024- Pavel Plotnitskii. All rights reserved.
/// This file is part of the \b stencil project.
///
/// \b stencil is free software: you may redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// The stencil project is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with the \b stencil project. If not, see <http://www.gnu.org/licenses/>.
///
/// @author Pavel Plotnitskii
/// @file stencil/parser.h
/// @brief A customized user options parser.
///
/// Parses user options from command line or from a text file.
///
#include <stdbool.h>
#include <stdlib.h>


/// @brief Describes the type of a user option
typedef enum __option_kind {
  BOOL     = 0,
  INT         ,
  FLOAT       ,
  STRING      ,
  VECT_INT    ,
  VECT_FLOAT  ,
  VECT_STRING
} option_kind;

/// @brief Typedef @ref _option for convenience
typedef struct __option option;

/// @brief A struct that describes a single user option
struct __option {
  option_kind kind;
  char    short_name;
  char   *long_name;
  char   *default_value;
  char   *value;
  char   *message;
  option *next;
};

/// @brief A linked list that contains all the user options
typedef struct __parser {
  char exe_name[128];
  char* description;
  option* options;
} parser;

/// @brief Creates a user option
/// @param kind the type of the option
/// @param short_name is a single character to reference the option
/// @param long_name a long name that characterizes the option
/// @param default_value is the value assigned to the option once created
/// @param message a help message that describes the option
/// @return a pointer to the option
option* option_create(const option_kind kind,
                      const char short_name,
                      const char *long_name,
                      const char *default_value,
                      const char *message);

/// @brief Deletes a user option
/// @param o pointer to the option to be deleted
void option_delete(option *o);

/// @brief Prints informations about a user option
/// @param o pointer to the option to be printed
void option_print(option *o);

/// @brief Checks if a user option is of @ref BOOL kind
/// @param o pointer to the option to be checked
bool option_boolean(option *o);

/// @brief Deletes a user option
/// @param o pointer to the option to be deleted
void option_delete(option *o);

/// @brief Prints informations about a user option
/// @param o pointer to the option to be printed
void option_print(option *o);

/// @brief Checks if a user option is of a @ref VECT_INT,
/// of a @ref VECT_STRING or of a @ref VECT_FLOAT kind
/// @param o pointer to the option to be checked
bool option_vect(option *o);

/// @brief Creates a parser
/// @param desc is a description of the executable
/// @return a pointer to the created parser
parser* parser_create(const char* desc);

/// @brief Deletes a parser
/// @param p pointer to the parser to be deleted
void parser_delete(parser* p);

/// @brief Prints all available options to the user
/// @param p pointer to the parser to be printed
void parser_usage(parser* p);

/// @brief Finds an option by its long name
/// @param p pointer to the parser that may contain the option
/// @param name the long name of the option
/// @return a pointer to he option if found or NULL
option* parser_find_option_by_long_name(parser* p, const char* name);

/// @brief Finds an option by its short name
/// @param p pointer to the parser that may contain the option
/// @param name the short name of the option
/// @return a pointer to he option if found or NULL
option* parser_find_option_by_short_name(parser* p, char name);

/// @brief Parses the command line entries
/// @param p pointer to a parser that will contain the parsed options
/// @param argc the number of command line entries
/// @param argv an array that contains the command line entries
/// @return the index of the last parsed entry
int  parser_parse(parser* p, int argc, char* argv[]);

/// @brief Looks for options in a text file.
/// The format of the text entries should be key=value
/// @param p pointer to a parser that will contain the parsed options
/// @param file_name the name of the text file to be parsed
/// @return the index of the last parsed entry
int  parser_parse_from_file(parser* p, char* file_name);

/// @brief Creates and adds a new option to an existing parser
/// @param p a pointer to the parser
/// @param kind the type of the option
/// @param short_name is a single character to reference the option
/// @param long_name a long name that characterizes the option
/// @param default_value is the default value assigned to the option once created
/// @param message a help message that describes the option
void parser_put(parser* p, const option_kind kind,
                const char short_name, const char *long_name,
                const char *default_value, const char *message);

/// @brief Checks if the user has asked help
/// @param p a pointer to the parser
bool  parser_help(parser* p);

/// @brief Detects the value of a @ref BOOL option
/// @param p a pointer to the parser
/// @param expression a string that contains the long name or
///  the short name of the option
/// @return the value of the option
bool  parser_get_bool(parser* p, const char* expression);

/// @brief Detects the value of a @ref INT option
/// @param p a pointer to the parser
/// @param expression a string that contains the long name or
///  the short name of the option
/// @return the value of the option
int   parser_get_int(parser* p, const char* expression);

/// @brief Detects the value of a @ref FLOAT option
/// @param p a pointer to the parser
/// @param expression a string that contains the long name or
///  the short name of the option
/// @return the value of the option
float parser_get_float(parser* p, const char* expression);

/// @brief Detects the value of a @ref STRING option
/// @param p a pointer to the parser
/// @param expression a string that contains the long name or
///  the short name of the option
/// @return the value of the option
char* parser_get_string(parser* p, const char* expression);

/// @brief Detects the values of a @ref VECT_INT option
/// @param p a pointer to the parser
/// @param expression a string that contains the long name or
///  the short name of the option
/// @param v a pointer to an array that will contains the values of the option
/// @param n the number of the values of the option
void  parser_get_vect_int(parser* p,
                          const char* expression, int* v, unsigned int n);

/// @brief Detects the values of a @ref VECT_INT option (size_t)
/// @param p a pointer to the parser
/// @param expression a string that contains the long name or
///  the short name of the option
/// @param v a pointer to an array that will contains the values of the option
/// @param n the number of the values of the option
void  parser_get_vect_size_t(parser* p,
                             const char* expression, size_t* v, unsigned int n);

/// @brief Detects the values of a @ref VECT_FLOAT option
/// @param p a pointer to the parser
/// @param expression a string that contains the long name or
///  the short name of the option
/// @param vf a pointer to an array that will contains the values of the option
/// @param n the number of the values of the option
void  parser_get_vect_float(parser* p,
                            const char* expression, float* vf, unsigned int n);

/// @brief Detects the values of a @ref VECT_STRING option
/// @param p a pointer to the parser
/// @param expression a string that contains the long name or
///  the short name of the option
/// @param vs a pointer to an array that will contains the values of the option
/// @param n the number of the values of the option
void  parser_get_vect_string(parser* p,
                             const char* expression, char** vs, unsigned int n);

/// @brief Sets the possible options for \b stencil
///
/// run <b> bin/stencil -h </b> or <b> bin/stencil --help </b>
/// to see the possible options
#define PARSER_BOOTSTRAP(p)                                             \
parser_put(p, BOOL,'v', "verbose", "false", "enable verbose mode");             \
parser_put(p, BOOL,'c', "cpu", "false", "run serial code on the CPU");          \
parser_put(p, VECT_INT,'l', "local", "16,4", "set the GPU block dimensions");       \
parser_put(p, INT,'d', "device", "0", "select the GPU device");                \
parser_put(p, BOOL,'o', "one", "false", "use only one GPU kernel");             \
parser_put(p, STRING,0, "gpu_options", " ", "set GPU build options");             \
parser_put(p, STRING,0, "in", "NONE", "velocity file");                           \
parser_put(p, STRING,0, "dir", "NONE", "output directory");                       \
parser_put(p, INT,0, "n1", "100", "velocity x dimension");                     \
parser_put(p, INT,0, "n2", "100", "velocity y dimension");                     \
parser_put(p, INT,0, "n3", "100", "velocity z dimension");                     \
parser_put(p, INT,0, "drcv", "10", "receivers step");                          \
parser_put(p, INT,0, "dshot", "5", "shot step");                               \
parser_put(p, INT,0, "dcdp", "10", "space delta step");                        \
parser_put(p, INT,0, "dline", "10", "space delta step");                       \
parser_put(p, INT,0, "ddepth", "10", "space delta step");                      \
parser_put(p, INT,0, "src_depth", "0", "source depth");                        \
parser_put(p, INT,0, "rcv_depth", "5", "receivers depth");                     \
parser_put(p, FLOAT,0, "dx", "10", "space delta step");                          \
parser_put(p, FLOAT,0, "dy", "10", "space delta step");                          \
parser_put(p, FLOAT,0, "dz", "10", "space delta step");                          \
parser_put(p, INT,'i', "iter", "1000", "simulation time step number");         \
parser_put(p, FLOAT,0, "dt", "0.001", "simulation time sampling in sec");         \
parser_put(p, FLOAT,0, "cfl", "0.8", "CFL percentage");                          \
parser_put(p, FLOAT,0, "fmax", "25.", "source max frequency");                   \
parser_put(p, FLOAT,0, "vmin", "1500.", "min velocity");                         \
parser_put(p, FLOAT,0, "vmax", "4500.", "max velocity");                         \
parser_put(p, INT,0, "nbsnap", "-1", "snapshot frequency");                    \
parser_put(p, BOOL,0, "check", "false", "check the GPU results");               \
parser_put(p, FLOAT,'e', "epsilon", "1.e-2", "the margin of floating point errors");\
parser_put(p, INT,0, "first", "-1", "first shot");                             \
parser_put(p, INT,0, "last",  "-1", "last shot");                              \
parser_put(p, INT,0, "tb_thread_group_size", "1","tb_thread_group_size");   \
parser_put(p, INT,0, "tb_nb_thread_groups", "1","tb_nb_thread_groups");  \
parser_put(p, INT,0, "tb_th_x", "1", "tb_th_x");                               \
parser_put(p, INT,0, "tb_th_y", "1", "tb_th_y");                               \
parser_put(p, INT,0, "tb_th_z", "1", "tb_th_z");                               \
parser_put(p, INT,0, "tb_t_dim", "1","tb_t_dim");                              \
parser_put(p, INT,0, "tb_num_wf", "1","tb_num_wf");                             \
parser_put(p, INT,0, "mode", "2","mode (1 MEM 2 I/O, 0 FUSE)");                     \
parser_put(p, INT,0, "fwd_steps", "1","tb param, save wavefield every fwd_steps diamonds");   \
parser_put(p, STRING,0,"tb_affinity", "NONE", "Affinity setup file");              \
parser_put(p, INT,0, "cbx", "10", "SB cache blocking in x");                               \
parser_put(p, INT,0, "cby", "22", "SB cache blocking in y");                               \
parser_put(p, INT,0, "cbz", "9999", "SB cache blocking in z");                               \
parser_put(p, INT,0, "rec_sismos", "1", "1 if record seismograms, 0 if not");                               \
parser_put(p, INT,0,"order", "1","solve acoustic wave equation of 1st or 2nd order, choose int 1 or 2");
#endif //  __STENCIL_PARSER_H_
