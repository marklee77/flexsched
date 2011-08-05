#!/bin/sh
(
    echo "set xlabel \"coefficient of variation\""
    echo "set ylabel \"minimum yield\""
    echo "set key bottom outside horizontal width +1" 
    #echo "set term png size 1024,768"
    echo "set term postscript eps color"
    for S in {1..9}; do
        echo "set title \"Minyield vs. Actual Coefficient of Variation for 64 hosts, 250 tasks, slack = 0.${S}\""
        echo "set output \"minyield-scatterplot-vs-acov-hosts-64-tasks-250-slack-0.$S.eps\""
        echo "plot for [i=0:1] \"scatterplot-0.$S.txt\" index i using 2:4 with points title columnheader(1)"
    done
) | gnuplot
