#!/bin/bash
#SBATCH --account=ucgd-kp
#SBATCH --partition=ucgd-kp
#SBATCH -o %j-%N.out
#SBATCH -e %j-%N.err
set -exo pipefail -o nounset
awk '$1 !~ /GL|NC|hs/ {print $0}' /scratch/ucgd/serial/quinlan_lab/data/u1021864/LCR-hs37d5.bed > poo.bed
dir=($(ls $DATA/wg*.bed*; ls $DATA/hg*; ls $DATA/LCR*))
i=$1
bedtools intersect -a $SCRATCH/domaincoordsdnds.bed -b <(awk '{if ($4<=0.1) print $0}' ${dir[$i]} | perl -pe 's/^chr(.*)$/$1/g') -f 0.2 -sorted > $SCRATCH/$(echo ${dir[$i]} | perl -pe 's/.*\/(.*)\..*/$1/g')'.txt'
#run with for loop, that's what i is for

./bedtools intersect -a $SCRATCH/domaincoordsdnds.bed -b $DATA/hginterrupted.bed $DATA/hgmicro.bed $DATA/hgrepeats.bed $DATA/hgsegmental.bed $DATA/hgselfchain.bed $DATA/hgsimple.bed $DATA/LCR-hs37d5.bed $SCRATCH/wgEncodeCrgMapabilityAlign100merscrubbed.bed $SCRATCH/wgEncodeCrgMapabilityAlign24merscrubbed.bed $SCRATCH/wgEncodeCrgMapabilityAlign36merscrubbed.bed $SCRATCH/wgEncodeCrgMapabilityAlign40merscrubbed.bed $SCRATCH/wgEncodeCrgMapabilityAlign50merscrubbed.bed $SCRATCH/wgEncodeCrgMapabilityAlign75merscrubbed.bed $SCRATCH/wgEncodeDacMapabilityConsensusExcludablescrubbed.bed $SCRATCH/wgEncodeDukeMapabilityRegionsExcludablescrubbed.bed $SCRATCH/wgEncodeDukeMapabilityUniqueness20bpscrubbed.bed $SCRATCH/wgEncodeDukeMapabilityUniqueness35bpscrubbed.bed -v -sorted > $SCRATCH/norepeats-maps.bed