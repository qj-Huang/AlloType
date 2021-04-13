# AlloType
Methods to predict allosteric regulation type.

AlloType -- A Coarsd-Grained Model to predict allosteric regulation type

### Installation Guide
First check that you have CMake version 3.1.0 or later. Check that you have compiled Boost and Eigen3.

For Ubuntu, you can install essential components by executing the following command:

```sh
sudo apt-get install cmake build-essential libboost-dev libboost-system-dev libboost-filesystem-dev libeigen3-dev
```

AlloType can be compiled from source code using CMake. Enter the working directory that contains this `README.md` file, then executing these commands:

```sh
mkdir build
cd build
cmake ..
make
sudo make install
```

AlloType will be installed into `/usr/bin` by default. You can specify another installation location by adding `--DCMAKE_INSTALL_PREFIX=/path/to/AlloType/root` to the end of cmake command.

There is one simple example in the `examples/` directory that go over the basic way to use AlloType.

For detailed information about the underlying mechanisms, please see [PAPER].

AlloType was written by Qiaojing Huang and Pengbo Song in the [Liu Lab](https://www.chem.pku.edu.cn/liuzhirong/) and [Lai Lab](http://mdl.ipc.pku.edu.cn/mdlweb/home-cn.php).

Please address all questions to chemhqj@pku.edu.cn 
