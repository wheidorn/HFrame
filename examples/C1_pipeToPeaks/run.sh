#/bin/bash

# in case of error, you should add the path like (depending on your OS):
#export DYLD_LIBRARY_PATH=$PWD/../bin:$DYLD_LIBRARY_PATH


make clean
make
if [[ $? != 0 ]]; then
  break
fi

mkdir -p ../../plot
while read file stavename temp side; do 

  if [[ "${file}" == "" ]]; then continue; fi
  if [[ ${file:0:1} == "#" ]]; then continue; fi

  
  echo 'Name:'$stavename'' > config
  echo 'Side:'$side'' >> config
  echo 'Tliquid:'$temp'' >> config
  echo running $stavename

  ./dana -name $stavename
  mkdir -p ../../plot/$stavename
  mv *.png *.pdf ../../plot/$stavename

done < ../../input.txt
