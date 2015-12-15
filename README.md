# Parallel Particle Mesh applied to N-body simulations

<p align="center"
![Evolution](https://github.com/sheroze1123/ppm/blob/master/particles.gif)
</p>

## Building ##
Our particle simulators depend on the FFTW discrete Fourier transform library.
To download, build, and install the library, simply run the following command
in the root `ppm` directory.

```
./scripts/fftw-install.sh
```

## Visualization ##
This repository contains two ways to visualize particle simulations:

1. You can generate an **mp4 file**. mp4 files are convenient for sharing
   visualizations with others and they let you step through a visualization
   frame by frame. The disadvantage is that mp4 files can take up quite a bit
   of space, and they cannot be generated in real-time (more on this later).
2. You can view a simulation in a **GUI**. Viewing a simulation in a GUI is
   nice because you don't need to make large mp4 files, and the GUI can animate
   simulations in real-time. The disadvantages are that the GUI is finicky to
   get working on all platforms, and you can't stop or slow down the animation.

Additionally, there are two sources used to generate a visualization.

1. **Files**. A particle simulation can write the positions of the particles
   into a text file that can be parsed by the visualizers. The advantage of the
   text files is that they can be used by either the mp4 or GUI visualizers.
   The disadvantage is that they get pretty large on big simulations.
2. **Network**. Alternatively, the simulators can stream data over the network
   to the GUI visualizer. The advantage of this is that the simulation can run
   for arbitrarily long without consuming a large amount of memory, and the
   simulations can be visualized in real-time! The downside is that the
   network can be slow, and mp4 files cannot be generated from streaming
   data.

In this section, we describe how to install all the dependencies you need to
get the visualizers working. Then, we discuss how to use the visualizers.

### Installation ###
Our two visualization scripts---[`scripts/mp4.py`](scripts/mp4.py) and
[`scripts/gui.py`](scripts/gui.py)---depend on a couple of packages and python
modules. Perform the following steps on your local machine to get them
installed. First, download a couple of packages needed for the GUI

```bash
sudo apt-get install tk-dev python-imaging-tk   # on linux
echo "TODO: not sure yet how to do this on OSX" # on OSX
```

Next, [install `pip`](https://pip.pypa.io/en/stable/) if you don't already have
it installed.  Then, install `virtualenv` if you don't already have it
installed (i.e. `pip install virtualenv`). Next, set up a virtual environment
and install the modules needed to run the visualizer.

```bash
virtualenv .                    # create a virtual environment
source bin/activate             # set up the virtual environment
pip install -r requirements.txt # install the necessary modules
```

Note that `requirements.txt` might be a bit too strict on requirements. If it's
giving you trouble, you can also try the following.

```bash
virtualenv .
source bin/activate
pip install numpy moviepy pillow
```

`numpy` and `moviepy` are needed by `scripts/mp4.py`; `numpy` and `pillow` are
needed by `scripts/gui.py`.

### Running the Scripts ###
Both scripts have a set of common flags. For example, both contain a `-p` flag
to set the number of pixels in the animation, a `-f` flag to set the FPS, and
`-m` flag to set the amount of mass needed to make a pixel completely white.
`mp4.py` and `gui.py` can both produce visualizations from a text file, and
`gui.py` can produce visualizations from over the network. Here are some
examples.

```bash
# create an mp4 (`particles.txt`)
python scripts/mp4.py particles.txt

# create an mp4 (`output.mp4`)
python scripts/mp4.py -o output.mp4 particles.txt

# create a 10 FPS mp4
python scripts/mp4.py -f 10 particles.txt

# create a 300x300 mp4
python scripts/mp4.py -p 300 particles.txt

# use the GUI with a file
python scripts/gui.py file particles.txt

# all the common flags still apply
python scripts/gui.py -f 10 -p 300 file particles.txt

# use the GUI from over the network
python scripts/gui.py net localhost 8000
```

TODO: explain how to get the gui working over the network.
