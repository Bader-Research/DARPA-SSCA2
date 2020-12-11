#ifndef _ALG_CREATE_INPUT_H
#define _ALG_CREATE_INPUT_H

#include "simple.h"

void create_input_nas(int n, int *x);
void create_input_nas_simple(int n, int *x, THREADED);

void create_input_random(int myN, int *x);
void create_input_random_simple(int myN, int *x, THREADED);

#endif
