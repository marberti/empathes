/****************************************************************************
 * Copyright (C) 2020-2022  Marco Bertini                                   *
 *                                                                          *
 * This file is part of Empathes.                                           *
 *                                                                          *
 * Empathes is free software: you can redistribute it and/or modify         *
 * it under the terms of the GNU General Public License as published by     *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 *                                                                          *
 * Empathes is distributed in the hope that it will be useful,              *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 * GNU General Public License for more details.                             *
 *                                                                          *
 * You should have received a copy of the GNU General Public License        *
 * along with Empathes.  If not, see <https://www.gnu.org/licenses/>.       *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/*** Function Declarations **************************************************/

void help(char *prog_name);
void open_file(const char *fname, FILE **fstream);
void close_file(const char *fname, FILE **fstream);
int  search_string(FILE *fstream, const char searched_str[],
        char returned_str[], size_t returned_str_len);
void get_best_iteration(FILE *fstream, int *iteration_n,
        double *highest_norm, int flag_spin);

/*** Main *******************************************************************/

int
main(int argc, char *argv[])
{
    char *fname;
    FILE *fstream;
    int str_n;
    int iteration_n;
    double highest_norm;
    const char searched_str[] = "PES Mode";
    const int flag_spin = 1;

    if (argc != 2) {
        help(argv[0]);
        exit(EXIT_FAILURE);
    } else if ((strcmp(argv[1],"-h")==0)||(strcmp(argv[1],"--help")==0)) {
        help(argv[0]);
        exit(EXIT_SUCCESS);
    }
    fname = argv[1];
    open_file(fname,&fstream);
    str_n = search_string(fstream,searched_str,NULL,0);
    if (str_n == 0) {
        printf("String \"%s\" not founded\n",searched_str);
        close_file(fname,&fstream);
        exit(EXIT_FAILURE);
    }
    printf("String \"%s\" founded on line %d\n",searched_str,str_n);
    get_best_iteration(fstream,&iteration_n,&highest_norm,flag_spin);
    printf("Best iteration     : %d\n",iteration_n);
    printf("Highest force norm : %e\n",highest_norm);
    close_file(fname,&fstream);
    exit(EXIT_SUCCESS);
}

/*** Function Definitions ***************************************************/

void
help(char *prog_name)
{
    const char *msg = 
    "Get informations about the NEB iteration closest to convergence,\n"
    "namely the highest norm on the total force of the best iteration.\n"
    "This script is intended to be used on a failed Empathes calculation.\n";
    printf("Usage: %s <file>\n",prog_name);
    printf("%s",msg);
}

/****************************************************************************/

void
open_file(const char *fname, FILE **fstream)
{
    if ((*fstream = fopen(fname,"r")) == NULL) {
        printf("ERR Cannot open %s\n", fname);
        exit(EXIT_FAILURE);
    }
}

/****************************************************************************/

void
close_file(const char *fname, FILE **fstream)
{
    if (fclose(*fstream) != 0) {
        printf("ERR Cannot close %s\n", fname);
        exit(EXIT_FAILURE);
    }
}

/****************************************************************************/

int
search_string(FILE *fstream, const char searched_str[],
        char returned_str[], size_t returned_str_len)
{

/****************************************************************************
 * This function search for the substring searched_str in each line of the  *
 * file associated with the file descriptor fstream. If a match is found,   *
 * the subroutine returns an integer with the number of the line where it   *
 * occurred, otherwise, 0 is returned.                                      *
 * If returned_str is not NULL, it will be used to store the first          *
 * returned_str_len characters of the line that generated a match, if any.  *
 ****************************************************************************/

    const int BUFF_LEN = 200;
    const int SEARCHED_LEN = strlen(searched_str);
    int line;
    int ch;
    int compared;
    char buff[BUFF_LEN+1];
    char original_buff[BUFF_LEN+1];
    
    buff[BUFF_LEN] = '\0';
    for (line = 1;; ++line) {
        for (int i = 0;; ++i) {
            ch = fgetc(fstream);
            if (ch == '\n') {
                if (i < BUFF_LEN) buff[i] = '\0';
                break;
            } else if (ch == EOF) {
                if (returned_str != NULL) {
                    returned_str = NULL;
                }
                return(0);
            } else if (i < BUFF_LEN) {
                buff[i] = (char) ch;
            }
        }
        strcpy(original_buff,buff);
        while (true) {
            if (strlen(buff) < SEARCHED_LEN) break;
            for (compared = 0; compared < SEARCHED_LEN; ++compared) {
                if (buff[compared] != searched_str[compared]) break;
            }
            if (compared == SEARCHED_LEN) {
                if (returned_str != NULL) {
                    strncpy(returned_str,original_buff,returned_str_len-1);
                    returned_str[returned_str_len-1] = '\0';
                }
                return(line);
            }
            for (int i = 0; i < strlen(buff); ++i) {
                buff[i] = buff[i+1];
            }
        }
    }
}

/****************************************************************************/

void
get_best_iteration(FILE *fstream, int *iteration_n, double *highest_norm,
    int flag_spin)
{

/****************************************************************************
 * This function reads all the iterations in the output file of Empathes    *
 * associated with the file descriptor fstream, and finds out which one is  *
 * the closest to convergence.                                              *
 * For each iteration, the highest norm on the total force (one for image)  *
 * is taken, being the best iteration the one whith the smallest,           *
 * highest norm. The best iteration is returned in iteration_n, while the   *
 * corresponding highest norm is returned in highest_norm.                  *
 * The flag_spin parameter states if the spin column is present (1) or      *
 * absent (0) in the empathes output file.                                  *
 ****************************************************************************/

    const char start_iteration[] = "Iteration";
    const char table_header[] = "Energy";
    int line;
    double norm;
    double max_norm;
    int image_n;
    int readed;

    *iteration_n = 0;
    *highest_norm = 0.0;
    for (int i = 1;; ++i) {
        line = search_string(fstream,start_iteration,NULL,0);
        if (line == 0) break;
        line = search_string(fstream,table_header,NULL,0);
        if (line == 0) break;
        if (i == 1) {
            for (int j = 0;; ++j) {
                if (flag_spin) {
                    readed = fscanf(fstream,"%d %*f %*s %lf %*s",&image_n,&norm);
                }
                else {
                    readed = fscanf(fstream,"%d %*f %lf %*s",&image_n,&norm);
                }
                if (readed < 2) break;
                if (j == 0) {
                    max_norm = norm;
                } else {
                    if (norm > max_norm) max_norm = norm;
                }
            }
            *iteration_n = i;
            *highest_norm = max_norm;
        } else {
            for (int j = 0; j < image_n; ++j) {
                if (flag_spin) {
                    readed = fscanf(fstream,"%*d %*f %*s %lf %*s",&norm);
                }
                else {
                    readed = fscanf(fstream,"%*d %*f %lf %*s",&norm);
                }
                if (readed < 1) {
                    printf("Error while reading iteration %d\n",i);
                    exit(EXIT_FAILURE);
                }
                if (j == 0) {
                    max_norm = norm;
                } else {
                    if (norm > max_norm) max_norm = norm;
                }
            }
            if (max_norm < *highest_norm) {
                *iteration_n = i;
                *highest_norm = max_norm;
            }
        }
    }
}
