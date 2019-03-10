#!/bin/bash

if [ "$1" != "" ]; then
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

  for pop in 10 100 1000
  do
    for file in $1/*_${pop}_*.dat
    do
      cat $file >> "pop_$pop"
    done
  done

  for mnm in 3 7 9
  do
    for file in $1/*_${mnm}_*.dat
    do
      cat $file >> "mnm_$mnm"
    done
  done

  for states in 6 18 24
  do
    for file in $1/*_${states}.dat
    do
      cat $file >> "states_$states"
    done
  done

else
  echo "Must provide name for directory"
fi





