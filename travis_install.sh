#!/bin/sh
set -e

# Install Eigen3 from source (for specified $EIGEN version)
if [ ! -d "$HOME/eigen/$EIGEN/include/eigen3/Eigen" ]; then
    wget https://bitbucket.org/eigen/eigen/get/$EIGEN.tar.bz2 -O /tmp/eigen.tar.bz2
    mkdir eigen3-src && tar -xvjf /tmp/eigen.tar.bz2 -C eigen3-src --strip-components 1
    mkdir -p $HOME/eigen/$EIGEN
    cd eigen3-src && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$HOME/eigen/$EIGEN .. && make install && cd ../..
else
    echo "Using cached Eigen3 $EIGEN installation";
fi
