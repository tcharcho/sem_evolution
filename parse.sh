#!/bin/bash

if [ "$1" != "" ]; then

  # ------------- Extract fitness data ------------- #
  for file in $1/*.sem
  do
    # echo $file
    filename="$(cut -d'/' -f2 <<<"$file")"
    # echo $filename
    newFile="$(cut -d'.' -f1 <<<"$filename")"
    # echo $newFile

    grep fit $file | sort -n > $1/$newFile.dat

    sed -i '' 's/ -fitness//g' $1/$newFile.dat
  done


  # ------------- group data based on hyperparameter type and value ------------- #

  # Population
  for pop in 5 2000 5000 #10 100 1000
  do
    rm "$1/pop_$pop"
    for file in $1/*_${pop}_*.dat
    do
      cat $file >> "$1/pop_$pop"
    done
  done

  # MNM
  for mnm in 3 7 9
  do
    rm "$1/mnm_$mnm"
    for file in $1/*_${mnm}_*.dat
    do
      cat $file >> "$1/mnm_$mnm"
    done
  done

  # Number of states
  for states in 3 50 70 #6 18 24
  do
    rm "$1/states_$states"
    for file in $1/*_${states}.dat
    do
      cat $file >> "$1/states_$states"
    done
  done

else
  echo "Must provide name for directory"
fi





