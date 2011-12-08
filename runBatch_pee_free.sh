#!/bin/bash
######### Is this datasets are pbpb or pp?
datasets="/afs/cern.ch/user/m/miheejo/scratch0/cms442_Jpsi/src/B2Jpsi/PEE/datasets"
output="fracfree"
executable=$1
errrange=0.01-0.2

rapbins=(0.0-2.4 0.0-1.2 1.2-1.6 1.6-2.4)
centbins=(0-10 10-20 20-30 30-40 40-50 50-100)
pt0012=(6.5-10.0 10.0-30.0 6.5-30.0)
pt1216=(5.5-10.0 6.5-10.0 10.0-30.0 5.5-30.0 6.5-30.0)
pt1624=(3.0-10.0 6.5-10.0 10.0-30.0 3.0-30.0 6.5-30.0)
pt0024=(0.0-30.0 6.5-10.0 10.0-30.0 6.5-30.0)


########## Function for progream running
function program {
  mc1="/afs/cern.ch/user/m/miheejo/scratch0/cms442_Jpsi/src/B2Jpsi/Jpsi_Histos_MC_NP.root"    #Non-prompt MC
  mc2="/afs/cern.ch/user/m/miheejo/scratch0/cms442_Jpsi/src/B2Jpsi/Jpsi_Histos_MC_PR.root"    #Prompt MC
  storage=/castor/cern.ch/user/m/miheejo/JpsiStudy/2011/70ub_result/
  executable=$1
  input=$2
  output=$3
  rap=$4
  cent=$5
  pt=$6
  err=$7
  ismb=$8
  file="$input/Data2011_cent$cent.root";
  work=$output"_rap"$rap"_cent"$cent"_pT"$pt;
  echo "Processing: "$work
  script="$executable -f $file -m $mc1 $mc2 -c -b -u -p $pt -y $rap -t $cent -e $err -l 1.0-2.0 -x $ismb -z 1 >& $work.log;";
  awk -v p=$(pwd) -v p2="$script" -v p3=$work -v p4=$storage '{gsub("_pwd_",p); gsub("&",p2); gsub("_output_file_",p3); gsub("_storage_",p4); print;}' < mjob.csh > $work.csh;
  bsub -R "pool>10000" -q 1nd -J $work < $work.csh
}


######### fracfree for minbias fitting
#for rap in ${rapbins[@]}
#do
#  cent=0-100  #Minbias fitting goes first
#  if [ "$rap" == "0.0-2.4" ]; then
#    for pt in ${pt0024[@]}
#    do
#      program $executable $datasets $output $rap $cent $pt $errrange 0
#    done
#  elif [ "$rap" == "0.0-1.2" ]; then
#    for pt in ${pt0012[@]}
#    do
#      program $executable $datasets $output $rap $cent $pt $errrange 0
#    done
#  elif [ "$rap" == "1.2-1.6" ]; then
#    for pt in ${pt1216[@]}
#    do
#      program $executable $datasets $output $rap $cent $pt $errrange 0
#    done
#  elif [ "$rap" == "1.6-2.4" ]; then
#    for pt in ${pt1624[@]}
#    do
#      program $executable $datasets $output $rap $cent $pt $errrange 0
#    done
#  fi
#done

######### fracfree fitting
if [ $1 == "Fit2DDatapp" ]; then
  exit 0
fi

for xopt in 1      # option for "x" : isMB(0), isMBPT(1), !isMBPT(2)
do
  for cent in ${centbins[@]}; do
    for rap in ${rapbins[@]}; do
      if [ "$rap" == "0.0-2.4" ]; then
        for pt in ${pt0024[@]}; do
          program $executable $datasets $output $rap $cent $pt $errrange $xopt
        done
      elif [ "$rap" == "0.0-1.2" ]; then
        for pt in ${pt0012[@]}; do
          program $executable $datasets $output $rap $cent $pt $errrange $xopt
        done
      elif [ "$rap" == "1.2-1.6" ]; then
        for pt in ${pt1216[@]}; do
          program $executable $datasets $output $rap $cent $pt $errrange $xopt
        done
      elif [ "$rap" == "1.6-2.4" ]; then
        for pt in ${pt1624[@]}; do
          program $executable $datasets $output $rap $cent $pt $errrange $xopt
        done
      fi
    done
  done

done
