from matplotlib import animation
from matplotlib import pyplot as plt
import numpy as np
from random import randint

def black(n, m):
    """
    Returns a n by m RGBA numpy array where each pixel is black with 100%
    opacity.
    """
    a = np.zeros((n, m, 4))
    a[:,:,3] += 1
    return a

def render(points, a):
    """
    Given a list of (x, y) coordinates (`points`) and a numpy array (`a`)
    decrease the opacity of a[x, y] for each (x, y) in `points`. This has the
    effect of whitening the position of each particle. If multiple particles
    overlap, the position becomes whiter. It is a precondition that all (x, y)
    coordinates are within `a`.
    """
    dalpha = 0.5
    for (x, y) in points:
        a[x,y,3] = min(a[x,y,3] - dalpha, 1)

def random_velocity():
    v = 2
    return (randint(-v, v), randint(-v, v))

def main():
    fps = 30                # frames per second
    seconds = 3             # mp4 length in seconds
    frames = fps * seconds  # number of mp4 frames
    n, m = 250, 250         # height (n) and width (n) of mp4 in pixels
    nparticles = n * m / 10 # number of particles

    fig = plt.figure()
    im = plt.imshow(black(n, m), interpolation="none")
    points = [(randint(0, m-1), randint(0, n-1)) for _ in range(nparticles)]
    vs = [random_velocity() for _ in range(nparticles)]

    def init():
        a = black(n, m)
        render(points, a)
        im.set_array(a)
        return [im]

    def animate(i):
        for (i, (x, y)) in enumerate(points):
            (dx, dy) = vs[i]
            points[i] = ((x + dx) % m, (y + dy) % n)
        a = black(n, m)
        render(points, a)
        im.set_array(a)
        return [im]

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=frames, blit=True)
    anim.save('particles.mp4', fps=fps)

if __name__ == "__main__":
    main()
