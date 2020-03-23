# Medit [![Build Status](https://travis-ci.org/ISCDtoolbox/Medit.svg?branch=master)](https://travis-ci.org/ISCDtoolbox/Medit)
OpenGL-based scientific visualization software

Medit was developped to visualize numerical simulation results on unstructured meshes in two and three dimensions. Scalar, vector and tensor fields can be easily associated and displayed with meshes.

#### Installation

In a terminal, clone this repository:

   ` git clone https://github.com/ISCDtoolbox/Medit.git `

   navigate to the downloaded directory:

   ` cd Medit `

   then create build directory and compile the project using cmake
   ```
   mkdir build
   cd build
   cmake ..
   make
   make install
   ```

#### Usage

* Simply run :
    `medit yourFile.mesh`

* Thanks to C. Dobrzynski, there is an [inline HTML documentation](https://www.ljll.math.upmc.fr/frey/logiciels/Docmedit.dir/index.html) available in french.

* There is also a [technical report](https://www.ljll.math.upmc.fr/frey/publications/RT-0253.pdf) describing its main features.

#### Remark :

* The experimental branch `feature/multi-sols` allows the visualisation of a solution field stored in a .sol file that contains multiple solution fields.
