# Parallel Particle Mesh applied to N-body simulations

## Visualization ##
First, build and run the serial implementation of the code:

```bash
make && ./serial
```

This will generate a file `particles.txt`. Run the visualizer with
`particles.txt` as input:

```bash
python visualizer.py particles.txt
```

This will generate `particles.mp4`. Running `visualizer.py` on totient might be
a bit slow. You may want to copy `particles.txt` to your host machine and run
the visualizer there.
