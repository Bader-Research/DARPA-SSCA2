#!/bin/sh
mpirun -np 4 -machinefile $PWD/machines.alpha.B -nolocal $PWD/fastest -m $PWD/machines.alpha.B -- $*
