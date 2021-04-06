//  Copyright (C) 2020-2021     Marco Bertini
//  
//  This file is part of neb.x.
//  
//  neb.x is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//  
//  neb.x is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with neb.x.  If not, see <https://www.gnu.org/licenses/>.

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

