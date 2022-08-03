#!/bin/bash
for i in plot*sc
do
  sed "s|data/Elise/Mistra_Elise_data/Elise/Mistra_Elise_|data/Elise/Mistra_Elise_|g" $i > tmp.sc
  sed "s|_2_unfinished_run||g" tmp.sc > $i
#  mv tmp.sc $i
  chmod u+x $i
done
