#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// function declarations

extern void f_error(const char *msg);
void c_chdir(char *dir);

// function definitions

void c_chdir(char *dir) {

  const int err_msg_len = 300;
  char err_msg[err_msg_len];

  //printf(" c_chdir: Changing directory to \"%s\"\n",dir);

  if (chdir(dir) != 0) {
    strcpy(err_msg,"c_chdir: cannot change to \"");
    strcat(err_msg,dir);
    strcat(err_msg,"\" directory");
    err_msg[err_msg_len-1] = '\0';
    f_error(err_msg);
  }

}

