#!/bin/bash

job=$1
test=0
if [[ $job < 0 ]]; then
    job=$((-job))
    test=1
fi
inpath=root://cmseos.fnal.gov//store/user/tucker/TrackMoverV1

fns=(
JetHT2016.root
JetHT2016B3.root
JetHT2016BCD.root
JetHT2016BCDEF.root
JetHT2016BthruG.root
JetHT2016C.root
JetHT2016D.root
JetHT2016E.root
JetHT2016EF.root
JetHT2016F.root
JetHT2016G.root
JetHT2016GH.root
JetHT2016H.root
JetHT2016H2.root
JetHT2016H3.root
qcdht1000.root
qcdht1000_hip1p0_mit.root
qcdht1500.root
qcdht1500_hip1p0_mit.root
)

nls=( 2 3 )
nbs=( 0 1 2 )

nfns=${#fns[@]}
nnls=${#nls[@]}
nnbs=${#nbs[@]}
njobs=$((nfns * nnls * nnbs))
echo inpath is $inpath
echo \#fns $nfns \#nls $nnls \#nbs $nnbs max jobs $njobs

nmax=$((nfns * nnls * nnbs))
echo job $job, nmax is $nmax
if [[ $job -ge $nmax ]]; then
    echo problem
    exit 1
fi

ifn=$((job % nfns))
fn=${fns[$ifn]}
iz=$((job / nfns))
inl=$((iz / nnbs))
inb=$((iz % nnbs))
nl=${nls[$inl]}
nb=${nbs[$inb]}
z=$nl$nb
outfn=${z}_$fn
echo fn $fn iz $iz inl $inl nl $nl inb $inb nb $nb z $z outfn $outfn

if [[ $test == 1 ]]; then
    echo test only, possibly modifying hists.jdl and hists_finish.sh
    sed -i -e "s/Queue.*/Queue ${njobs}/" hists.jdl
    sed -i -e "s/^njobstot=.*/njobstot=${njobs}/" hists_finish.sh
    exit 1
fi

export SCRAM_ARCH=slc6_amd64_gcc530
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW CMSSW_8_0_25 2>&1
cd CMSSW_8_0_25/src
eval $(scram ru -sh)
cd ../..

./hists.exe $inpath/$z/$fn $outfn $nl $nb 2>&1
