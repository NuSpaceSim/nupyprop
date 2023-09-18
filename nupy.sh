#!/bin/bash

set -e

export PATH
mkdir nupyenv
tar -xzf nupyprop.tar.gz -C nupyenv 
. nupyenv/bin/activate

ENERGY=$1
#ANGLE=$2
STATS=$2
PN_MODEL=$3
HTC=$4
JOB=$5

#nupyprop -e $ENERGY -a $ANGLE -s $STATS -htc $HTC -job $JOB
#nupyprop -e $ENERGY -a $ANGLE -s $STATS -p $PN_MODEL -htc $HTC -job $JOB
nupyprop -e $ENERGY -s $STATS -p $PN_MODEL -htc $HTC -job $JOB
