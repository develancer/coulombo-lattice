#!/bin/sh
cd `dirname $0`
echo "-- EXPECTED --"
cat <<TXT
 1  1  1  1   0.017643339  0.000000000
 1  1  2  2  -0.000026347  0.000040578
 1  2  2  1   0.017413492  0.000000000
 2  2  2  2   0.017726404  0.000000000
TXT
echo "-- ACTUAL --"
../coulombo --dielectric=12.4 --spin --step=10 \
            --integrals=1111,1122,1221,2222 \
            e1-D.arma e1-U.arma h2-D.arma h2-U.arma
