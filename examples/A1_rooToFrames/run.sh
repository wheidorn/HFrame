#/bin/bash

# in case of error, you should add the path like (depending on your OS):
#export DYLD_LIBRARY_PATH=$PWD/../bin:$DYLD_LIBRARY_PATH

make clean
make
if [[ $? != 0 ]]; then
  break
fi

mkdir -p ../../roo/combined

while read file stavename temp side; do 

  if [[ "${file}" == "" ]]; then continue; fi
  if [[ ${file:0:1} == "#" ]]; then continue; fi

  
  echo 'Name:'$stavename'' > config
  echo 'Side:'$side'' >> config
  echo 'Tliquid:'$temp'' >> config
  echo running $stavename

  ./dana -name $file
  mv out.txt ../../roo/combined/${file}.txt #Moved both of these to the main folder
  ls -l ../../roo/combined/${file}.*

done < ../../input.txt
