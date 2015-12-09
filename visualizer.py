from matplotlib import animation
from matplotlib import pyplot as plt
import numpy as np
import random

def black():
    """
    Returns a 100x100 RGBA numpy array where each pixel is black with 100%
    opacity.
    """
    a = np.zeros((100, 100, 4))
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
    dalpha = 1
    for (x, y) in points:
        a[x,y,3] = min(a[x,y,3] - dalpha, 1)

def main():
    n, m = 100, 100
    fig = plt.figure()
    im = plt.imshow(black(), interpolation="none")

    npoints = 50
    points = [(random.randint(0, n-1), random.randint(0, n-1)) for _ in range(npoints)]
    v = 2

    def init():
        a = black()
        render(points, a)
        im.set_array(a)
        return [im]

    def animate(i):
        points[:] = [((x + v) % m, (y + v) % n) for (x, y) in points]
        a = black()
        render(points, a)
        im.set_array(a)
        return [im]

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=256, interval=20, blit=True)
    anim.save('basic_animation.mp4', fps=10)

if __name__ == "__main__":
    main()
