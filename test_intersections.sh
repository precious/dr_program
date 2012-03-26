#!/bin/bash

OUTPUT_FILE="test_intersections.result"
NUMBER_OF_TESTS=100
NUMBER_OF_INTERVALS=20


echo -n > $OUTPUT_FILE
for i in `seq 1 $NUMBER_OF_TESTS` ; do
	echo -n "$i "
	./program -t 30000 out.bnd >> $OUTPUT_FILE
done
echo

DATA_FILE="`./generate_gnuplot_data.py \"$OUTPUT_FILE\" $NUMBER_OF_INTERVALS`"
IMAGE_FILE="freq.png"
gnuplot -e "set term png; set output \"$IMAGE_FILE\"; plot \"$DATA_FILE\" using 1:2 with lines;"
xdg-open "$IMAGE_FILE"
