# Parallel Particle Mesh applied to N-body simulations

## Visualization ##
Our visualization script, [`scripts/visualizer.py`](scripts/visualizer.py),
depends on a couple of python modules that are not installed on the totient
cluster. Perform the following steps on your local machine. First, [install
`pip`](https://pip.pypa.io/en/stable/) if you don't already have it installed.
Then, install `virtualenv` if you don't already have it installed (i.e. `pip
install virtualenv`). Next, set up a virtual environment and install the
modules needed to run the visualizer.

```bash
virtualenv .                    # create a virtual environment
source bin/activate             # set up the virtual environment
pip install -r requirements.txt # install the necessary modules
```

Note that `requirements.txt` might be a bit too strict on requirements. If it's
giving you trouble, you can also try the following.

```bash
virtualenv .              # create a virtual environment
source bin/activate       # set up the virtual environment
pip install numpy moviepy # install the necessary modules
```

Next, get a simulation text file, (e.g. `particles.txt`) from the totient
cluster and run it through the visualizer.

```bash
python scripts/visualizer.py particles.txt
```

This will generate `particles.mp4`.
