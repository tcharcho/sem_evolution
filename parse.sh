#!/bin/bash
for file in full_run_1/*.sem
do
  # echo $file
  filename="$(cut -d'/' -f2 <<<"$file")"
  # echo $filename
  newFile="$(cut -d'.' -f1 <<<"$filename")"
  # echo $newFile

  grep fit $file | sort -n > full_run_1/$newFile.dat
done




