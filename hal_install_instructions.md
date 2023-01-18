# hal & HALPER installation instructions

Instructions to install [hal](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md) and [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing). 

These have been tested on Mac OS Monterey 12.6. These instructions should work on Unix systems as well.

1. Install Homebrew (requires root access/`sudo`).

*Not necessary on Unix systems.*

Why: so that you can install `wget`.

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
This will take ~10 minutes. It might seem to get stuck at "Receiving objects" for a few minutes; this is normal, it should finish eventually.

After this script runs, make sure to follow the `Next steps` to finish the install:
```
==> Next steps:
- Run these three commands in your terminal to add Homebrew to your PATH:
    echo '# Set PATH, MANPATH, etc., for Homebrew.' >> /Users/<username>/.zprofile
    echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> /Users/<username>/.zprofile
    eval "$(/opt/homebrew/bin/brew shellenv)"
```

2. Install `wget`.

*Not necessary on Unix systems.*

Why: so that you can download the `hdf5` installer.

```
brew install wget
```

3. Download `hdf5`.

```
cd ~/Downloads/
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz --no-check-certificate
```

4. Clone [hal](https://github.com/ComparativeGenomicsToolkit/hal/blob/master/README.md) and [HALPER](https://github.com/pfenninglab/halLiftover-postprocessing) repos.

These instructions will assume that you are cloning your repos into a directory called `~/repos/`. If you clone your repos into a different path, please update any paths that begin with `~/repos/` accordingly.
```
mkdir ~/repos
cd ~/repos
git clone https://github.com/ComparativeGenomicsToolkit/hal.git
git clone https://github.com/pfenninglab/halLiftover-postprocessing.git
```

5. Install Anaconda (`conda`).

If `conda` is already installed on your system, you can skip this step.

```
# We will download the install file, run it, and then delete it in a temporary directory
cd /tmp

# Please choose the correct installer for your OS and architecture.
# To determine your architecture on Mac, use "About This Mac"
# all installers listed at https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links
# Linux 64-bit: Miniconda3-latest-Linux-x86_64.sh
# macOS Apple M1 64-bit: Miniconda3-latest-MacOSX-arm64.sh
# macOS Intel x86 64-bit: Miniconda3-latest-MacOSX-x86_64.sh

# Replace Miniconda3-latest-Linux-x86_64.sh with the file you chose in the previous step
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the Miniconda installer. You will need to answer a few questions. The default options should be fine.
bash Miniconda3-latest-Linux-x86_64.sh

# Delete the installer and return to the repos dir
rm Miniconda3-latest-Linux-x86_64.sh
cd ~/repos
```

6. Create `hal` conda env.

```
conda create -n hal python=3 numpy matplotlib
```

7. Configure `hdf5`.

```
cd ~/Downloads
tar xzf hdf5-1.10.1.tar.gz
cd hdf5-1.10.1
mkdir ~/repos/hal/hdf5
./configure --enable-cxx --prefix ~/repos/hal/hdf5
```

The output should look something like this:
```
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
```

8. Build `hdf5` inside the `hal` repo.

```
cd ~/Downloads/hdf5-1.10.1
make && make install
```

9. Clone and build [sonLib](https://github.com/ComparativeGenomicsToolkit/sonLib):

```
cd ~/repos
git clone https://github.com/ComparativeGenomicsToolkit/sonLib.git
pushd sonLib && make && popd
```

The output should look something like this:
```
~/repos/sonLib ~/repos
[...]
cc -Iinc -Iimpl -fsigned-char -I../externalTools/quicktree_1.1/include/ -I ../externalTools/cutest -fPIC -std=c99 -fsigned-char -O3 -g -Wall --pedantic -funroll-loops -DNDEBUG   -o ../../sonLib/bin/sonLib_fastaCTest tests/fastaCTest.c ../../sonLib/lib/sonLib.a    -lz -lm 
cp sonLib_daemonize.py ./bin/sonLib_daemonize.py
chmod +x ./bin/sonLib_daemonize.py
~/repos
```
It is ok if you get several warnings.

10. Set environment variables to build `hal`.

For these commands, please replace `[repos dir]` with the full path to your `~/repos` directory, e.g. `/Users/<username>/repos`. If you use `~/repos`, it is not guaranteed to work.

You do not need to add these lines to your `~/.bashrc` or `~/.bash_profile`. These environment settings are for building `hal` only.

```
export PATH=[repos dir]/hal/hdf5/bin:${PATH}
export h5prefix=-prefix=[repos dir]/hal/hdf5
```

11. Build `hal`:

```
cd ~/repos/hal
make 
```

12. Setup environment variables to run `hal` and `HALPER`:

Add the following to your `~/.bash_profile`. Please replace `[repos dir]` with the full path to your `~/repos` directory, as in step 9.

```
# hal and HALPER
export PATH=[repos dir]/hal/bin:${PATH}
export PYTHONPATH=[repos dir]/halLiftover-postprocessing:${PYTHONPATH}
```

13. Test that you can run `hal`:

```
source ~/.bash_profile
halLiftover
```

Expected output:
```
Too few (required positional) arguments

halLiftover v2.2: Map BED or PSL genome interval coordinates between two genomes.

USAGE:
halLiftover [Options] <halFile> <srcGenome> <srcBed> <tgtGenome> <tgtBed>
[...]
```

14. Test that you can run `HALPER`:
```
conda activate hal
python -m orthologFind
```

Expected output:
```
usage: orthologFind.py [-h] [-max_len MAX_LEN] [-max_frac MAX_FRAC]
                       [-protect_dist PROTECT_DIST] [-min_len MIN_LEN]
                       [-min_frac MIN_FRAC] -qFile QFILE -tFile TFILE -sFile
                       SFILE -oFile OFILE [-mult_keepone] [-narrowPeak]
                       [-noHist] [-keepChrPrefix KEEPCHRPREFIX]
orthologFind.py: error: the following arguments are required: -qFile, -tFile, -sFile, -oFile
```
