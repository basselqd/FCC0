# FCC Nanosphere Simulation Scripts

The above scripts are used to perform ompi simulations of FCC array of "atomistic" maghemite nanoparticles using the MagLua scripting environment.

MagLua is series of extensions to the Lua scripting language (lua-5.1.5 from https://www.lua.org/ftp/lua-5.1.5.tar.gz). MagLua was developed by Jason Mercer to provide a set of computational tools used to facilitate the numerical simulation of magnetic systems using a variety of different techniques. The resultant code is very efficient as most of the subroutines are written in C++ and can utilise the Message Passing Interface (MPI) when appropriate, to exploit the benefits of parallelism.

In the course of its development MagLua has undergone a number of revisions that have extended its capabilities. All versions of MagLua are available online:

https://github.com/jasonimercer/maglua

svn versions have been converted to git tags and are available as the following:

https://github.com/jasonimercer/maglua/tree/r365

The above scripts require MagLua r-307. This may be downloaded from

https://github.com/jasonimercer/maglua/releases/tag/r307

The MagLua extensions are well documented. The file `maglua.html` contains more information about MagLua and the documentation specific to Revision-307.

_________________________________________________________________________________

The above files comprise the following:

`MKT.lua` will generate the matrices (`8x8x8_fcc.lua`) used to calculate  the energy and and the effective fields for an 8x8x8 FCC lattice of point dipoles with periodic boundary conditions.

`A.lua` is the main file where the 8x8x8 dipole array is constructed from the magnetic moments of each of the individual nanoparticles and the dipole fields acting on them calculated. The value of the dipole fields is then communicated to each of the nanoparticle on the FCC lattice using MPI.

`DenseSphere.lua`contains the scripts for an initial run (ie. not from a previously saved point) to build objects and functions for the nanoparticles.

`DenseSphereL.lua` is similar to  `DenseSphere.lua` but is used when the system is loaded from saved point.

`SetupSphericalParticle.lua` contain the functions required to build the nanoparticle crystal with the specified exchange and anisotropy parameters.

_________________________________________________________________________________

With the Lua scripting language and the MagLua extensions installed, a simulation with some default parameters is initiated by entering a command as:
````sh
$ mpirun -np 512 maglua A.lua
````

Initiating a simulation from a saved point at temperature $T so that it runs to final time $t enter the command

````sh
$ mpirun -np 512 maglua A.lua load $T again $t
````

The files used to store the data used to recreate the saved point have names like "SSS$T.$rank.dat"


