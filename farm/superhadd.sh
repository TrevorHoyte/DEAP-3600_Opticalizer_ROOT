#!/bin/bash
#SBATCH --account=rrg-deap
#SBATCH --time=04:00:00 
#SBATCH --mem-per-cpu=8G 
#SBATCH --cpus-per-task=1 
PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH
DEAPHOME=/project/6004969
PATH=$PATH:$DEAPHOME/software/bin:
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DEAPHOME/software/lib
MANPATH=$MANPATH:$DEAPHOME/software/share/man
export DEAPHOME PATH LD_LIBRARY_PATH MANPATH MACHTYPE OSTYPE
PATH=$PATH:$HOME/.local/bin:$HOME/bin
source $DEAPHOME/software/ratcage/env.sh

cd /scratch/trevorh/rat6_ambe/ntp/ 
hadd -f ~/data/ambe_mc_rat6.root *ntp.root
echo done

#this file adds all the ntp files together into one big file. 
#edit lines 16 to match  your in the diretcory with all the ntp 
# edit line 17 to match the desired location and name of the combined file 