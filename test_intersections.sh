#!/bin/bash

OUTPUT_FILE="test_intersections.result"
NUMBER=100

echo -n > $OUTPUT_FILE
for i in `seq 1 $NUMBER` ; do
	echo -n "$i "
	./program -t 3000000 out.bnd >> $OUTPUT_FILE
done
echo
