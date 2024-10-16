/// @file src/parser.c
/// @brief This file contains the implementation of the option parser.
///
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stencil/parser.h>
#define MSG(fmt,...)                    \
  do {                                  \
    fprintf(stdout, "[STENCIL MSG]:");  \
    fprintf(stdout, fmt, ##__VA_ARGS__);\
    fprintf(stdout, "\n");              \
  } while(0)

option *option_create(const option_kind kind,
                      const char short_name,
                      const char *long_name,
                      const char *default_value,
                      const char *message) {
  option *o        = (option*) malloc(sizeof(option));
  o->kind          = kind;
  o->short_name    = short_name;
  o->long_name     = (char*) malloc(strlen(long_name)+1);
  o->default_value = (char*) malloc(128);
  o->value         = (char*) malloc(128);
  o->message       = (char*) malloc(strlen(message)+1);
  o->next          = NULL;
  strcpy(o->long_name,     long_name);
  strcpy(o->default_value, default_value);
  strcpy(o->value,         default_value);
  strcpy(o->message,       message);
  return o;
}

void option_delete(option *o) {
  free(o->long_name);
  free(o->default_value);
  free(o->value);
  free(o->message);
  free(o);
}

void option_print(option* o) {
  fprintf(stdout, "\t%c%c%c --%-12s %s   %s.\n",
    o->short_name == 0 ? ' ' : '-',
    o->short_name == 0 ? ' ' : o->short_name,
    o->short_name == 0 ? ' ' : ',',
    o->long_name,
    option_boolean(o) ? "         ":
    option_vect(o) ? "\033[4mV0,V1,...\033[0m" :
    "\033[4mVALUE\033[0m    ",
    o->message);
}

bool option_boolean(option* o) { return o->kind == BOOL; }
bool option_vect(option* o) {
  return
  o->kind == VECT_INT   ||
  o->kind == VECT_FLOAT ||
  o->kind == VECT_STRING;
}

parser* parser_create(const char* description) {
  parser *p      = (parser*) malloc(sizeof(parser));
  p->description = (char*) malloc(strlen(description)+1);
  p->options     = option_create(BOOL, 'h', "help", "false",
                                 "print this help message");
  strcpy(p->description, description);
  parser_put(p, STRING, 'f', "file",  "", "parse parameters from file");
  return p;
}

void parser_delete(parser* p) {
  option *n, *o = p->options;
  do {
    n = o->next;
    option_delete(o);
    o = n;
  } while(o != NULL);
  free(p->description);
  free(p);
}

void parser_usage(parser *p) {
  option *o = p->options;
  fprintf(stdout, "\nNAME\n");
  fprintf(stdout, "\t%s - %s", p->exe_name, p->description);
  fprintf(stdout, "\nSYNTAX\n");
  fprintf(stdout, "\t%s [options]\n", p->exe_name);
  fprintf(stdout, "\nOPTIONS\n");
  while(o->next != NULL) option_print(o = o->next);
  fprintf(stdout, "\n");
}

option* parser_find_option_by_long_name(parser* p, const char* name) {
  option *o = p->options;
  while(o != NULL) {
    if (strcmp(o->long_name, name) == 0) break;
    o = o->next;
  }
  return o;
}

option* parser_find_option_by_short_name(parser* p, char name) {
  option *o = p->options;
  while(o != NULL) {
    if ((o->short_name) && (o->short_name == name)) break;
    o = o->next;
  }
  return o;
}

char *get_basename(char *path) {
  char *base = strrchr(path, '/');
  return base ? base+1 : path;
}

int parser_parse(parser *p, int argc, char* argv[]) {
  option *o;
  unsigned int s, i, param_size;
  const char *param;
  char *long_name, short_name;
//  MSG("long_name=%s,short_name=%s",long_name,short_name);
  strcpy(p->exe_name, get_basename(argv[0]));
  for (i = 1; i < argc; ++i) {
    param      = argv[i];
//    MSG("param=%s",param);
//    printf("param=%s\n",param);
    param_size = strlen(param);
    if (param[0] != '-') {
      if (((i+1) < argc) && (argv[i+1][0] != '-'))
        return i+1;
      continue;
    } else {
      if (param[1] == '-') {
        long_name = strndup(param+2, param_size-2);
        if ((o=parser_find_option_by_long_name(p, long_name)) == NULL) {
          fprintf(stderr,
            "[PARSER ERROR]: Option '%s' not recognized\n", long_name);
          exit(EXIT_FAILURE);
        }
        if (option_boolean(o)) {
          if((i+1) < argc && argv[i+1][0] != '-') {
            fprintf(stderr,
              "[PARSER ERROR]: Boolean '%s' option has to be unary\n",
              long_name);
            exit(EXIT_FAILURE);
          }
          strcpy(o->value, "true");
        } else {
          if((i+1) >= argc || strlen(argv[i+1]) == 0) {
            fprintf(stderr,
              "[PARSER ERROR]: '%s' option value missing or empty\n",
              long_name);
            exit(EXIT_FAILURE);
          }
          strcpy(o->value, argv[i+1]);
        }
        if (strcmp(long_name,"file") == 0) parser_parse_from_file(p, argv[i+1]);
        free(long_name); long_name = NULL;
      } else {
        for (s = 1; s < param_size; s++) {
          short_name = param[s];
          if ((o=parser_find_option_by_short_name(p, short_name)) == NULL) {
            fprintf(stderr, "[PARSER ERROR]: Option '%c' not recognized\n",
              short_name);
            exit(EXIT_FAILURE);
          }
          if (option_boolean(o)) {
            if((s == param_size-1) &&
             ((i+1) < argc)     &&
             (argv[i+1][0] != '-')) {
              fprintf(stderr,
                "[PARSER ERROR]: Boolean '%c' option has to be unary\n",
                short_name);
          }
          strcpy(o->value, "true");
        } else {
          if((s > 1) && (s < param_size-1)) {
            fprintf(stderr, "[PARSER ERROR]: Illegal options sequence\n");
            exit(EXIT_FAILURE);
          }
          if((i+1) >= argc || strlen(argv[i+1]) == 0) {
            fprintf(stderr,
              "[PARSER ERROR]: '%c' option value missing or empty\n",
              short_name);
            exit(EXIT_FAILURE);
          }
          strcpy(o->value, argv[i+1]);
        }
        if (short_name == 'f') parser_parse_from_file(p, argv[i+1]);
      }
    }
  }
}
if (parser_help(p)) {
  parser_usage(p);
  exit(EXIT_SUCCESS);
}
return i;
}


int parser_parse_from_file(parser *p, char* file_name) {
  FILE *f;
  option *o;
  char fline[128], *key, *value;
  if ((f = fopen(file_name, "r")) == 0) {
    fprintf(stderr, "[PARSER ERROR]: Couldn't open config file '%s'\n",
      file_name);
    exit(EXIT_FAILURE);
  }
  while (fgets(fline, 128, f) != NULL) {
    if (fline[0] == '#') {
      continue;
    } else {
      key   = strtok(fline, " =\n");
      value = strtok(0," =\n");
      if(strtok(0," =\n")) {
        fprintf(stderr,
          "[PARSER ERROR]: More then 2 entries per line are detected\n");
        exit(EXIT_FAILURE);
      }
      if(key) {
        if (strlen(key) == 1) {
                    // short option
          if(key[0] == 'f') {
            fprintf(stderr,
              "[PARSER ERROR]: Recursive file parsing is not supported\n");
            exit(EXIT_FAILURE);
          }
          if ((o=parser_find_option_by_short_name(p, key[0])) == NULL) {
            fprintf(stderr, "[PARSER ERROR]: Option '%c' not recognized\n",
              key[0]);
            exit(EXIT_FAILURE);
          }
          if (option_boolean(o)) {
            if(value) {
              fprintf(stderr,
                "[PARSER ERROR]: Boolean '%c' option has to be unary\n",
                key[0]);
            }
            strcpy(o->value, "true");
          } else {
            if(value == NULL) {
              fprintf(stderr,
                "[PARSER ERROR]: '%c' option value missing or empty\n",
                key[0]);
              exit(EXIT_FAILURE);
            }
            strcpy(o->value, value);
          }
        } else {
                    // long option
          if(strcmp(key, "file") == 0) {
            fprintf(stderr,
              "[PARSER ERROR]: Recursive file parsing is not supported\n");
            exit(EXIT_FAILURE);
          }
          if ((o=parser_find_option_by_long_name(p, key)) == NULL) {
            fprintf(stderr,
              "[PARSER ERROR]: Option '%s' not recognized\n", key);
            exit(EXIT_FAILURE);
          }
          if (option_boolean(o)) {
            if(value) {
              fprintf(stderr,
                "[PARSER ERROR]: Boolean '%s' option has to be unary\n",
                key);
              exit(EXIT_FAILURE);
            }
            strcpy(o->value, "true");
          } else {
            if(value == NULL) {
              fprintf(stderr,
                "[PARSER ERROR]: '%s' option value missing or empty\n",
                key);
              exit(EXIT_FAILURE);
            }
            strcpy(o->value, value);
          }
        }
      }
    }
  }
  fclose(f);
  if (parser_help(p)) {
    parser_usage(p);
    exit(EXIT_SUCCESS);
  }
  return 0;
}

void parser_put(parser* p,
  const option_kind kind,
  const char short_name,
  const char* long_name,
  const char* default_value,
  const char* message) {
  option *o = p->options;
//  MSG("long_name=%s",long_name);
  if((strlen(long_name) == 0) || (isspace(long_name[0]))) {
    fprintf(stderr, "[PARSER ERROR]: Blank Option detected\n");
    exit(EXIT_FAILURE);
  }
  if(parser_find_option_by_long_name(p, long_name) != NULL) {
    fprintf(stderr, "[PARSER ERROR]: Option '%s' already exists\n", long_name);
    exit(EXIT_FAILURE);
  }
  if(parser_find_option_by_short_name(p, short_name) != NULL) {
    fprintf(stderr, "[PARSER ERROR]: Option '%c' already exists\n", short_name);
    exit(EXIT_FAILURE);
  }
  while(o->next != NULL) o = o->next;
  o->next = option_create(kind, short_name,
    long_name, default_value, message);
  return;
}


bool parser_get_bool(parser *p, const char* param_name) {
  option *o;
  if (strlen(param_name) == 1) {
    if((o=parser_find_option_by_short_name(p , param_name[0])) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Boolean option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  } else {
    if((o=parser_find_option_by_long_name(p , param_name)) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Boolean option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  }
  return((strcmp(o->value, "true") == 0) || (strcmp(o->value, "1") == 0));
}

int parser_get_int(parser *p, const char* param_name) {
  option *o;
  if (strlen(param_name) == 1) {
    if((o=parser_find_option_by_short_name(p , param_name[0])) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Int option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  } else {
    if((o=parser_find_option_by_long_name(p , param_name)) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Int option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  }
  return atoi(o->value);
}

float parser_get_float(parser* p, const char* param_name) {
  option *o;
  if (strlen(param_name) == 1) {
    if((o=parser_find_option_by_short_name(p , param_name[0])) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Float option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  } else {
    if((o=parser_find_option_by_long_name(p , param_name)) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Float option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  }
  return (float)atof(o->value);
}

char* parser_get_string(parser *p, const char* param_name) {
  option *o;
  if (strlen(param_name) == 1) {
    if((o=parser_find_option_by_short_name(p , param_name[0])) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: String option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  } else {
    if((o=parser_find_option_by_long_name(p , param_name)) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: String option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  }
  return o->value;
}

void parser_get_vect_int(parser* p,
                         const char* param_name, int* vi, unsigned int n) {
  option *o;
  char *pch, temp[128];
  unsigned int i=0;
  if (strlen(param_name) == 1) {
    if((o=parser_find_option_by_short_name(p , param_name[0])) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Vect int option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  } else {
    if((o=parser_find_option_by_long_name(p , param_name)) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Vect int option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  }
  pch = strstr(o->value, ",");
  while (pch) {
    i++;
    pch = strstr(pch + 1, ",");
  }
  if( n != (i+1) ) {
    fprintf(
      stderr,
      "[PARSER ERROR]: Size of int vector should be %d (actual = %d)\n",
      n, i+1);
    exit(EXIT_FAILURE);
  }
  strcpy(temp, o->value);
  pch = strtok(temp, ",");
  i = 0;
  while (pch != NULL) {
    vi[i++] = atoi(pch);
    pch = strtok(NULL, ",");
  }
}

void parser_get_vect_size_t(parser* p,
                            const char* param_name, size_t* vi, unsigned int n){
  option *o;
  char *pch, temp[128];
  unsigned int i=0;
  if (strlen(param_name) == 1) {
    if((o=parser_find_option_by_short_name(p , param_name[0])) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Vect int option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  } else {
    if((o=parser_find_option_by_long_name(p , param_name)) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Vect int option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  }
  pch = strstr(o->value, ",");
  while (pch) {
    i++;
    pch = strstr(pch + 1, ",");
  }
  if( n != (i+1) ) {
    fprintf(
      stderr,
      "[PARSER ERROR]: Size of int vector should be %d (actual = %d)\n",
      n, i+1);
    exit(EXIT_FAILURE);
  }
  strcpy(temp, o->value);
  pch = strtok(temp, ",");
  i = 0;
  while (pch != NULL) {
    vi[i++] = atol(pch);
    pch = strtok(NULL, ",");
  }
}

void parser_get_vect_float(parser *p,
                           const char* param_name, float* vf, unsigned int n) {
  option *o;
  char *pch, temp[128];
  unsigned int i=0;
  if (strlen(param_name) == 1) {
    if((o=parser_find_option_by_short_name(p , param_name[0])) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Vect float option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  } else {
    if((o=parser_find_option_by_long_name(p , param_name)) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Vect float option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  }
  pch = strstr(o->value, ",");
  while (pch) {
    i++;
    pch = strstr(pch + 1, ",");
  }
  if ( n != (i+1)) {
    fprintf(stderr, "[PARSER ERROR]: Size of float vector should be %d\n", n);
    exit(EXIT_FAILURE);
  }
  strcpy(temp, o->value);
  pch = strtok(temp, ",");
  i = 0;
  while (pch != NULL) {
    vf[i++] = (float)atof(pch);
    pch = strtok(NULL, ",");
  }
}

void parser_get_vect_string(parser *p,
                            const char* param_name, char** vs, unsigned int n) {
  option *o;
  char *pch, temp[128];
  unsigned int i=0;
  if (strlen(param_name) == 1) {
    if((o=parser_find_option_by_short_name(p , param_name[0])) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Vect string option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  } else {
    if((o=parser_find_option_by_long_name(p , param_name)) == NULL) {
      fprintf(stderr,
        "[PARSER ERROR]: Vect string option '%s' not recognized\n",
        param_name);
      exit(EXIT_FAILURE);
    }
  }
  pch = strstr(o->value, ",");
  while (pch) {
    i++;
    pch = strstr(pch + 1, ",");
  }
  if ( n != (i+1)) {
    fprintf(stderr, "[PARSER ERROR]: Size of string vector should be %d\n", n);
    exit(EXIT_FAILURE);
  }
  strcpy(temp, o->value);
  pch = strtok(temp, ",");
  i = 0;
  while (pch != NULL) {
    vs[i++] = strdup(pch);
    pch = strtok(NULL, ",");
  }
}

bool parser_help(parser *p) { return parser_get_bool(p, "help"); }

/// @}
