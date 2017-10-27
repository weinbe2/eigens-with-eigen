#! /bin/bash
# This script downloads Eigen and puts it in the directory specified by the first command line argument. If no argument is specified, it's downloaded to the current directory.

RELEVANT_DIR="./"

if [ "$#" -ge 1 ]; then
  RELEVANT_DIR="$1"
fi

cd $RELEVANT_DIR

if [ -d "Eigen" ]; then
  if [ -L "Eigen" ]; then
    echo "Eigen static link exists! Please delete it."
    exit
  else
    echo "Eigen directory exists! Please delete it."
    exit
  fi
fi

wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2

# When you download Eigen, it appends some arbitrary string to the end.
# We need to deal with.

TAR_NAME="3.3.4.tar.bz2"

mv $TAR_NAME eigen.tar.bz2
tar -xvjf eigen.tar.bz2

DIR_NAME=""
for f in eigen-*
do
  DIR_NAME=${f};
done

rm -rf Eigen
mv $DIR_NAME Eigen
rm eigen.tar.bz2


