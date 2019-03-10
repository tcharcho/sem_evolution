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
  done
else
  echo "Must provide name for directory"
fi





