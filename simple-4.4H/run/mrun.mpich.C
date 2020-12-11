#!/bin/sh
mpirun -np 2 -machinefile $PWD/machines.alpha.C -nolocal $PWD/fastest -m $PWD/machines.alpha.C -- $*

