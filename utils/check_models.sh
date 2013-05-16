IFS=$'\n'
DIR='models2'
for file in `ls $DIR`; do
    res="`./program -v -x $DIR/$file 2> /dev/null | grep polygons`"
    if [ $? == 0 ] ; then
        echo $file $res
    fi
done
