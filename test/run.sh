#/bin/bash

#export DYLD_LIBRARY_PATH=[HFRAME_DIR]/HFrames/bin:$DYLD_LIBRARY_PATH

gmake clean
gmake
if [[ $? != 0 ]]; then
  break
fi

mkdir -p roo/combined
mkdir -p roo/cropped

while read file stavename temp side; do 

  if [[ "${file}" == "" ]]; then continue; fi
  if [[ ${file:0:1} == "#" ]]; then continue; fi

  echo 'Name:'$stavename'' > config
  echo 'Side:'$side'' >> config
  echo 'Tliquid:'$temp'' >> config

  ./dana -name $file

done < input.txt

#
# the input.txt should be like
#
# <file name> <output name> <temperature setup> <L or J side>
#
# input assumed in roo/<file name>/seq_*.root
#longFlawRevNo4_sL_80lens_p50  StaveNo4_L_80lens_p50   50 L
