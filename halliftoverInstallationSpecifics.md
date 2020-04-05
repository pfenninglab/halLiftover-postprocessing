### Installation

For the most part, this follows the instructions at [ComparativeGenomicsToolkit/hal](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md). Some additional clarifications have been added.  These instructions have been tested using hdf5 versions 5-1.10.1 and 5-1.10.4; gcc versions 4.8.5, 4.9.2, and 5.3.0; and HAL Format API version 2.1.  These instructions were last updated on April 2, 2020 and are not being maintained.


### Downloading HAL Format API

From the parent directory of where you want HAL installed:
git clone git://github.com/glennhickey/hal.git


### Installing Dependencies

* Download hdf5:
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz --no-check-certificate

* Follow instructions for a 'Local install from source into DIR (do not need root password)'
	
	mkdir [directory with hal]/hal/hdf5

	tar xzf hdf5-1.10.1.tar.gz

	cd hdf5-1.10.1

	./configure --enable-cxx --prefix [directory with hal]/hal/hdf5

* output will end with something like this after installation:

				Features:
				---------
				                  Parallel HDF5: no
				             High-level library: yes
				                   Threadsafety: no
				            Default API mapping: v110
				 With deprecated public symbols: yes
				         I/O filters (external): deflate(zlib)
				                            MPE: no
				                     Direct VFD: no
				                        dmalloc: no
				 Packages w/ extra debug output: none
				                    API tracing: no
				           Using memory checker: no
				Memory allocation sanity checks: no
				            Metadata trace file: no
				         Function stack tracing: no
				      Strict file format checks: no
				   Optimization instrumentation: no

* Make hdf5 by running: 
		


		make && make install

		Output message will be something like this: 
				make[2]: Leaving directory `[directory with hal]/hal/hdf5-1.10.1/hl/examples'
				make[2]: Entering directory `[directory with hal]/hal/hdf5-1.10.1/hl/c++'
				make[3]: Entering directory `[directory with hal]/hal/hdf5-1.10.1/hl/c++/examples'
				../../../bin/mkdirs [directory with hal]/hal/hdf5/share/hdf5_examples/hl/c++
				+ /usr/bin/install -c ./ptExampleFL.cpp [directory with hal]/hal/hdf5/share/hdf5_examples/hl/c++/.
				+ /usr/bin/install -c run-hlc++-ex.sh [directory with hal]/hal/hdf5/share/hdf5_examples/hl/c++/.
				make[3]: Leaving directory `[directory with hal]/hal/hdf5-1.10.1/hl/c++/examples'
				make[2]: Leaving directory `[directory with hal]/hal/hdf5-1.10.1/hl/c++'
				make[1]: Leaving directory `[directory with hal]/hal/hdf5-1.10.1/hl'

				[[username]@[node-name] hdf5-1.10.1]$ cd [directory with hal]/hal/hdf5
				[[username]@[node-name] hdf5]$ ls
				bin  include  lib  share


### Prepare GCC 

Use gcc version >= 4.2.  If gcc has not been installed, install it.  If gcc has been installed as a module on a cluster, run:
	module load gcc[- gcc version]


### sonLib

From the same parent directory where you downloaded HAL:

	  git clone https://github.com/ComparativeGenomicsToolkit/sonLib.git
	  pushd sonLib && make && popd

NOTE: The version of zlib available in anaconda does not seem to be compatible with sonLib.


### Building HAL
* Update environment variables according to their instructions: 
		export PATH=[directory with hal]/hal/hdf5/bin:${PATH}
		export h5prefix=-prefix=[directory with hal]/hal/hdf5
* run the following (takes awhile):
		cd [directory with hal]/hal
		make 


### Running HAL

* add as given the following to your path:

		export PATH=[directory with hal]/hal/bin:${PATH}

		export PYTHONPATH=[directory with hal]/bin:${PYTHONPATH}


### ~/.bashrc additions that are useful for running hal:

```
export PATH=[directory with hal]/hdf5/bin:${PATH}
export h5prefix=-prefix=[directory with hal]/hal/hdf5
export sonLibRootPath=[directory with hal]/sonLib
export PATH=[directory with hal]/hal/bin:${PATH}
export PYTHONPATH=[directory with hal]:${PYTHONPATH}
export PYTHONPATH=[directory with hal]/bin:${PYTHONPATH}

```
