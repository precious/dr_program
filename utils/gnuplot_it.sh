#!/bin/sh

############## function for parsing command line arguments ##############
function __get_var_name {
    # arguments: OPTION, ARRAY
    # ARRAY -- each element is like '-f|--foo|...|{FOO}',
    # where -f, --foo, ... -- names of options
    #   FOO -- name of environment variable to be set
    # OPTION -- option like '-f' or '--foo'
    # function returns variable name for given option
    option="$1"
    shift
    should_stop=""
    # assumed that IFS is set to space
    for args in "$@"; do
        parts=(${args//|/ }) # here parts will be array like ("-f" "--foo" "FOO")
        for part in ${parts[@]:0:$[${#parts[@]} - 1]}; do
            if [ "$part" == "$option" ]; then
                # return last element without leftmost and rightmost characters
                echo ${parts[$[${#parts[@]} - 1]]:1:-1}
                should_stop="yes"
                break
            fi
        done
        if [ -n "$should_stop" ]; then
            break
        fi
    done
}

function echo_help {
    echo $1 | grep -Po '(?<=\[)[^]]+(?=\])' | sed 's/|{[^}]*}//' | sed "s/^/    /"
    echo "    -h|--help print help"
}

function opt2env {
    __DOCSTRING=$1
    shift # skip docstring
    __PREV_IFS="$IFS" # backup IFS
    IFS=$'\n' # set new IFS
    __PARAMS=`echo $__DOCSTRING | grep -Po '(?<=\[)[^]]+(?=\])'`
    req_arg=() # array of options that require arguments
    not_req_arg=() # array of options that don't require arguments
    for __P in $__PARAMS; do
        words=(${__P// /$'\n'}) # split param by spaces
        if (( ${#words[@]} > 0 )); then
            if (( ${#words[@]} > 1 )); then
                if echo ${words[1]} | grep -P '^<[^>]+>$' > /dev/null 2>1; then # option requires argument
                    req_arg+=("${words[0]}")
                else # option doesn't require argument
                    not_req_arg+=("${words[0]}")
                fi
            else # option doesn't require argument
                not_req_arg+=("${words[0]}")
            fi
        fi
    done
    
    IFS=' '
    RETVAL=0
    FREE_ARGUMENTS=()
    while [ $# -gt 0 ] ; do
        if [ "$1" == "-h" -o "$1" == "--help" ]; then
            echo_help "$__DOCSTRING"
            exit
        fi
        if [ "${1:0:1}" == "-" ]; then
            varname=`__get_var_name "$1" "${not_req_arg[@]}"`
            if [ -n "$varname" ]; then
                eval "$varname=1"
                shift
                continue
            fi
            varname=`__get_var_name "$1" "${req_arg[@]}"`
            if [ -n "$varname" ]; then
                eval "$varname=\"$2\""
                shift 2
                continue            
            fi
            # we encountered options that was not declared
            RETVAL=1 # set exit code to 1 and exit
            break
        else
            FREE_ARGUMENTS+=("$1")
            shift
        fi
    done
    IFS="$__PREV_IFS" # restore IFS
    # set exit code
    [ "$RETVAL" == "0" ] && true || false
}
#########################################################################

format=png
columns='1:2'
size='1300,700'

OPTSTRING="[-f|--format|{format} <format> -- output file format (default $format)]\
[-c|--columns|{columns} <columns> -- columns to plot (default $columns)]\
[-s|--size|{size} <size> -- size of output image (default $size)]"
opt2env "$OPTSTRING" "$@"
set -- "${FREE_ARGUMENTS[@]}"

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

