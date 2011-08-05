#!/bin/sh
(
    echo "set xlabel \"coefficient of variation\""
    echo "set ylabel \"average minimum yield\""
    echo "set xtics 0.1"
    echo "set key bottom outside horizontal width +1" 
    echo "set term png size 1024,768"
    for S in {1..9}; do
        echo "set title \"Average Minyield vs. Coefficient of Variation for 64 hosts, 250 tasks, slack = 0.${S}\""
        echo "set output \"minyield-vs-cov-hosts-64-tasks-250-slack-0.$S-hcpu.png\""
        echo "plot for [i=0:6] \"minyield-vs-cov-hosts-64-tasks-250-slack-0.$S-hcpu.txt\" index i using 1:4:5:6 with yerrorlines title columnheader(3)"
    done
    echo "set ylabel \"failure rate (%)\""
    echo "set yrange [0:100]"
    for S in {1..9}; do
        echo "set title \"Failure Rate vs. Coefficient of Variation for 64 hosts, 250 tasks, slack = 0.${S}\""
        echo "set output \"failrate-vs-cov-hosts-64-tasks-250-slack-0.$S-hcpu.png\""
        echo "plot for [i=0:6] \"minyield-vs-cov-hosts-64-tasks-250-slack-0.$S-hcpu.txt\" index i using 1:(100-\$3) with lines title columnheader(3)"
    done
) | gnuplot
