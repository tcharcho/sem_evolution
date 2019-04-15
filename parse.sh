#!/bin/bash

dirNames=("FF1_run_1" "FF2_run_1" "FF2_run_2" "production_run_1" "production_run_2" "production_run_3")

echo -e "\nExtracting marginal distrinution of fitness based on 3 hyperparameters:"
echo "- Population size"
echo "- Number of states"
echo -e "- Maximum number of mutations\n\n"

for dir in ${dirNames[@]}
do
  echo -e "Parsing data files in" $dir "...\n"

  pop=()
  mnm=()
  states=()

  for file in $dir/*.sem
  do
    # echo $file
    filename="$(cut -d'/' -f2 <<<"$file")"
    # echo $filename
    newFile="$(cut -d'.' -f1 <<<"$filename")"
    # echo $newFile

    IFS='_' read -ra ADDR <<< "$newFile"

    pop+=(${ADDR[1]})
    mnm+=(${ADDR[2]})
    states+=(${ADDR[3]})


  done

  pop=($(echo "${pop[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
  mnm=($(echo "${mnm[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
  states=($(echo "${states[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))


  # ------------- group data based on hyperparameter type and value ------------- #

  # Population
  for p in ${pop[@]}
  do
    rm "$dir/pop_$p.marg"
    for file in $dir/*_${p}_*.dat
    do
      cat $file >> "$dir/pop_$p.marg"
    done
  done

  # MNM
  for m in ${mnm[@]}
  do
    rm "$dir/mnm_$m.marg"
    for file in $dir/*_${m}_*.dat
    do
      cat $file >> "$dir/mnm_$m.marg"
    done
  done

  # Number of states
  for s in ${states[@]}
  do
    rm "$dir/states_$s.marg"
    for file in $dir/*_${s}.dat
    do
      cat $file >> "$dir/states_$s.marg"
    done
  done

done


