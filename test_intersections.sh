#!/bin/bash

NUMBER_OF_TESTS=700
NUMBER_OF_INTERVALS=10
NUMBER_OF_PARTICLES=30000
RADIUS=136
OUTPUT_FILE="${NUMBER_OF_TESTS}x${NUMBER_OF_PARTICLES}.result"

if [ ! -f "$OUTPUT_FILE" ] ; then
	echo -n > $OUTPUT_FILE
	for i in `seq 1 $NUMBER_OF_TESTS` ; do
		echo -n "$i "
		./program -t $NUMBER_OF_PARTICLES -r $RADIUS out.bnd >> "$OUTPUT_FILE"
	done
	echo
fi

DATA_FILE="`./generate_gnuplot_data.py \"$OUTPUT_FILE\" $NUMBER_OF_INTERVALS`"
IMAGE_FILE="${NUMBER_OF_TESTS}x${NUMBER_OF_PARTICLES}.png"
gnuplot -e "set term png; set output \"$IMAGE_FILE\"; plot \"$DATA_FILE\" using 1:2 with lines;"
rm -rf "$DATA_FILE"
xdg-open "$IMAGE_FILE"
