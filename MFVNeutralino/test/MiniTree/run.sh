listOfSamples="
background
mfv_neu_tau000100um_M0400_2017
mfv_neu_tau000100um_M0600_2017
mfv_neu_tau000100um_M0800_2017
mfv_neu_tau000100um_M1200_2017
mfv_neu_tau000100um_M1600_2017
mfv_stopdbardbar_tau000100um_M0400_2017
mfv_stopdbardbar_tau000100um_M0600_2017
mfv_stopdbardbar_tau000100um_M0800_2017
mfv_stopdbardbar_tau000100um_M1200_2017
mfv_stopdbardbar_tau000100um_M1600_2017
"

for sample in $listOfSamples
do
  for ntk in 3 4 5
  do
    python plot_triggers.py $sample $ntk
  done
done
