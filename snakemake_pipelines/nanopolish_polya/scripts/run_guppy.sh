#BSUB -J run_guppy
#BSUB -n 20
#BSUB -o %J.stdout
#BSUB -e %J.stderr
#BSUB -R span[hosts=1]
#BSUB -q gpu

export PATH=/work/bio-mowp/software/ont-guppy-6.0.0/ont-guppy/bin:$PATH

guppy_basecaller -i fast5 -s guppy_out -c rna_r9.4.1_70bps_hac.cfg --device "cuda:all:100%"
