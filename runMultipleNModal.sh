file=split_2d_modal.f90
nSTART=4
nOLD=$nSTART
for i in {4..7};
do

nNEW=$i
echo $nOLD $nNEW
echo "%s/polyOrder = $nOLD/polyOrder = $nNEW/g
w
q
" | ex $file
nOLD=$nNEW
make 2d_test
make clean

# Move output file to appropriate directory
mv dgnolim/spltMod2d_def_smth_cosbell.nc dgnolim/cnvg/n$nOLD/spltMod2d_def_smth_cosbell.nc
mv dgmfill/spltMod2d_def_smth_cosbell.nc dgmfill/cnvg/n$nOLD/spltMod2d_def_smth_cosbell.nc


done

echo "%s/polyOrder = $nOLD/polyOrder = $nSTART/g
w
q
" | ex $file

