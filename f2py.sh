#!/bin/bash

LDFLAGS="$LDFLAGS -shared" f2py -c --f90flags="-I../lib" ../lib/lib{traj,comp,varlist,ephrd,util,sys,navbase,math90}.a -m ssd ssd.f90   

wait

mv ssd*.so ssd.so

wait

cp ssd.so ../test/.

wait

echo 'all done'
exit 0


