#!/bin/sh

format=png
columns='1:2'
size='1300,700'

OPTIND=1
while getopts :hc:s:f: opt; do
    case "$opt" in
    h)
        echo "$0 [-c m:n] [-s 1300,700] [-f png] filename"
        exit 0
        ;;
    c)  columns="$OPTARG"
        ;;
    s)  size="$OPTARG"
        ;;
    f) format="$OPTARG"
        ;;

    esac
done
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

cat "$1" | grep -P '(\-?[0-9.]+\s+){5,}' > "${1}.filtered"
echo filtered

echo "set terminal $format
set output \"$1.$format\"
set terminal $format size $size
set format y \"%.12f\"
set format x \"%f\"
plot \"${1}.filtered\" using $columns with lines title \"заряд\"
quit
" | gnuplot

xdg-open "$1.$format"

