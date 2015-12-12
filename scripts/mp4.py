import argparse
import argparser
import math
import moviepy.editor as mpy
import numpy as np
import renderer

from parser import FileParser, NetParser
from random import randint
from renderer import Renderer

def main(args):
    parser = FileParser(args.filename)
    renderer = Renderer(args.pixels, args.mass, parser)

    fps = args.fps                    # frames per second
    frames = parser.n                 # number of mp4 frames
    duration = float(frames) / fps    # duration of mp4 in seconds
    nparticles = parser.N_p + 1       # number of particles

    clip = mpy.VideoClip(lambda _: renderer.next(), duration=duration)
    clip.write_videofile("particles.mp4", fps=fps)

def parse_args():
    parser = argparse.ArgumentParser()
    argparser.renderer_args(parser)
    argparser.file_args(parser)
    parser.add_argument(
        "-o", "--output",
        default="particles.mp4",
        help="output mp4 filename"
    )
    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())
