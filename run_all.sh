#!/bin/sh

FILES="$1/problem-*.txt"
for f in $FILES
do
  ./flexsched "$2" $f $f.result
done

