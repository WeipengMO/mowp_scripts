#BSUB -J run_guppy
#BSUB -n 8
#BSUB -o %J.stdout
#BSUB -e %J.stderr
#BSUB -R span[hosts=1]
#BSUB -gpu "num=2"


time guppy_basecaller \
  -i fast5 \
  -s guppy_out \
  -c dna_r9.4.1_450bps_hac.cfg \
  --recursive \
  --disable_pings \
  --qscore_filtering \
  --device "cuda:all:100%" 


