#include <stdio.h>
#include <string.h>
#include <error.h>
#include <stdbool.h>
#include <ctype.h>

void tolowerstr(char str[]);

int main(int argc, char *argv[])
{
  const int BUFF_LEN = 301;

  FILE *file_stream;
  char file_name[100];
  char subroutine_name[100];
  char buff[BUFF_LEN];
  char field[3][100];
  int  field_n;
  bool is_a_subroutine;
  bool subroutine_found;

  // cmd line argument section ----------------------------
  if (argc != 3) {
    error(1,0,"Required 2 arguments");
  }

  strcpy(file_name,argv[1]);
  strcpy(subroutine_name,argv[2]);
  tolowerstr(subroutine_name);

  // open file --------------------------------------------
  file_stream = fopen(file_name,"r");
  if (file_stream == NULL) {
    error(1,0,"Cannot open %s file",file_name);
  }

  // write subroutine to stdout ---------------------------
  is_a_subroutine  = false;
  subroutine_found = false;
  for (;;) {
    if (fgets(buff,BUFF_LEN,file_stream) == NULL) {
      break;
    }

    // read string from file
    field[0][0] = '\0';
    field[1][0] = '\0';
    field[2][0] = '\0';
    sscanf(buff,"%s %s %s",field[0],field[1],field[2]);
    tolowerstr(field[0]);
    tolowerstr(field[1]);
    tolowerstr(field[2]);

    // write subroutine body
    if (subroutine_found) {
      printf("%s",buff);
    }

    // search for subroutine start
    if (subroutine_found == false) {
      if ((strcmp(field[0],"subroutine") == 0) || (strcmp(field[0],"function") == 0)) {
        is_a_subroutine = true;
        field_n = 1;
      }
      else if (strcmp(field[1],"function") == 0) {
        is_a_subroutine = true;
        field_n = 2;
      }
      else {
        is_a_subroutine = false;
      }

      if (is_a_subroutine) {
        for (int i = 0; i < strlen(field[field_n]); i++) {
          if (field[field_n][i] == '(') {
            field[field_n][i] = '\0';
            break;
          }
        }

        if (strcmp(field[field_n],subroutine_name) == 0) {
          subroutine_found = true;
          printf("%s",buff);
        }
      }
    }
    // search for subroutine end
    else { // subroutine_found == true
      if ((strcmp(field[0],"end") == 0) &&
          ((strcmp(field[1],"subroutine") == 0) || (strcmp(field[1],"function") == 0))) {
        for (int i = 0; i < strlen(field[2]); i++) {
          if (field[field_n][i] == '(') {
            field[field_n][i] = '\0';
            break;
          }
        }

        if (strcmp(field[2],subroutine_name) == 0) {
          break;
        }
        else {
          error(1,0,"Bad subroutine ending");
        }
      }
    }
  }

  // close file -------------------------------------------
  if (fclose(file_stream) != 0) {
    error(1,0,"Cannot close %s file",file_name);
  }

  return 0;
}

//===================================================================

void tolowerstr(char str[])
{
  for (int i = 0; str[i] != '\0'; i++) {
    str[i] = tolower(str[i]);
  }
}
