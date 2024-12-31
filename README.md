Private Information Retrieval (PIR)
 Protocol faster Spiral
=====

***

 
The [faster Spiral](https://github.com/parsear/fspiral) is the implement of the paper [Faster Spiral: Low-Communication, High-Rate Private Information Retrieval]. This library will be continuously updated to support more features.

## Requirements
The runnig of this code requires a basic c++ toolchain, including

* g++ compiler (clang++12 or newer is optimal)
* make
* cmake


## Build and Run
We adapt [Intel HEXL](https://github.com/intel/hexl) library to implement the NTTs.
Therefore, please install it firstly.
We recommend that install it to your default user's program path `/usr/local`, and the suggeted version of them are showed in following table. What's more, using `clang++12` or newer compiler can achieve better performance.

| Libraries | Version | Website |
| ---- | ---- | ---- |
| HEXL  | v1.2.5 | https://github.com/intel/hexl |


We recommend using the following command to install HEXL. You can also refer to the website for more detailed information.

```
git clone https://github.com/intel/hexl.git -b v1.2.5 --single-branch 
cd hexl
mkdir build
cd build
cmake ..
make
sudo make install
```

After install [Intel HEXL](https://github.com/intel/hexl), you can build and run our PIR protocol by

```
cd fspiral
mkdir build
cd build
cmake ..
make
./tests/test-fspiral
```


An example output:
```
Database: 128 * 256 * R_p, (total 288 MB)

target (row, col) = (46, 229)
 preprocess costs 2553399 us.
keyGen Done.
 query generate costs 1959 us.
 online server response costs 1239194 us.
The err is [758268062, 2372643108, 66974689664262136, 2329731913, 676875373, 1281929247, 2376661815, 824912486, 66974689737105656, 1182736430, 1140187905, 1140267775, 1279705346, 1337510525, 1080414943, 961909084, ..., 66974687046202803]
 recover costs 275 us.

The decrypted value is [103937, 4728, 162395, 84485, 107484, 245835, 144626, 248980, 97776, 108081, 221164, 13151, 93770, 152874, 195204, 115399, ..., 159144]
The correct   value is [103937, 4728, 162395, 84485, 107484, 245835, 144626, 248980, 97776, 108081, 221164, 13151, 93770, 152874, 195204, 115399, ..., 159144]
```
