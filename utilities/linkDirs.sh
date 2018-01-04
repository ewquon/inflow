#!/bin/bash
#
# Utility to symlink time directories within constant/BCname/ for
#   periodic inflows
#
# written by Eliot Quon (eliot.quon@nrel.gov) -- 2017-07-21
#
# Note: Only handles integers for now

if [ -z "$1" ]; then
    echo "USAGE:"
    echo "  $0 [tEnd]"
    echo "  $0 [tStart] [tEnd]"
    exit
fi

tStart=0
tEnd="$1"

shift
if [ -n "$1" ]; then
    tStart="$tEnd"
    tEnd="$1"
fi

tStart0=
tEnd0=
tSkip=
for d in `ls -1vd *`; do
    if [ -d "$d" ]; then
        if [ -z "$tStart0" ]; then
            tStart0="$d"
        elif [ -z "$tSkip" ]; then
            tSkip=$((d-tStart0))
        else
            tEnd0="$d"
        fi
    fi
done
tLen=$((tEnd0-tStart0+tSkip))

echo "Start/end/delta/length time: $tStart0, $tEnd0, $tSkip, $tLen"

# DRY RUN FIRST
tgt=$tStart
src=$tStart
while [ "$tgt" -le "$tEnd" ]; do
    if [ ! -d "$tgt" ]; then
        if [ "$src" -ge "$tLen" ]; then
            src=$((src-tLen))
        fi
        echo "$tgt --> $src"
    fi
    src=$((src+tSkip))
    tgt=$((tgt+tSkip))
done
echo    '-----------------------------------------------'
read -p 'If this looks right, press enter to continue...'
echo    '-----------------------------------------------'

echo 'Making symlinks...'
tgt=$tStart
src=$tStart
while [ "$tgt" -le "$tEnd" ]; do
    if [ ! -d "$tgt" ]; then
        if [ "$src" -ge "$tLen" ]; then
            src=$((src-tLen))
        fi
        #ln -sv $src $tgt
        ln -s $src $tgt
    fi
    src=$((src+tSkip))
    tgt=$((tgt+tSkip))
done
echo 'Done.'
