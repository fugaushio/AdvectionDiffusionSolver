rm -r build
mkdir build
cd build
cmake -D compiler=intel \
      -D enable_GLOG=ON \
      ..

make && make install