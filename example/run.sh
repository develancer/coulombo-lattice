#!/bin/sh
cd `dirname $0`
../coulombo --atoms=X.3d --dielectric=11.4 --orbitals=10 --skip-lines=1 \
  e1.dat 2>/dev/null
cut -c 1-29 eeee.txt > eeee.cut
( diff - eeee.cut <<TXT
 1  1  1  1    0.090340071977
TXT
) && echo OK
rm -f eeee.*
