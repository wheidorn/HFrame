#/bin/bash


# in case of error, you should add the path like (depending on your OS):
#export DYLD_LIBRARY_PATH=$PWD/../bin:$DYLD_LIBRARY_PATH


gmake clean
gmake
if [[ $? != 0 ]]; then
  break
fi

mkdir -p plot
mkdir -p roo/Pipes
while read file stavename temp side; do 

  if [[ "${file}" == "" ]]; then continue; fi
  if [[ ${file:0:1} == "#" ]]; then continue; fi

  
  echo 'Name:'$stavename'' > config
  echo 'Side:'$side'' >> config
  echo 'Tliquid:'$temp'' >> config
  echo running $stavename

  ./dana -name $file
  mkdir -p plot/$stavename
  mv *.png plot/$stavename

  mv a1_pipes.root roo/Pipes/${stavename}.root
done < input.txt
