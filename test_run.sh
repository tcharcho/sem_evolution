#!/bin/sh

# popsize
for i in 5 500 2000 #10 100 1000
do
  # MNM
  for j in 3 7 9
  do
    # states
    for k in 3 10 30 #6 18 24
    do
      g++ -lm -O3 SemDistEvo.cpp sem.cpp stat.cpp -Dpopsize=$i -DMNM=$j -DSTATES=$k
      ./a.out $1
    done
  done
done
