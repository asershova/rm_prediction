#!/bin/bash
j=0
for i in $@; do
  
  tmux new-session -d -s my_session_$j "./run_specific_bash.sh $i"
  echo tmux new-session -d -s my_session_$j "./run_specific_bash.sh $i"
  j=$((j+1))
done

