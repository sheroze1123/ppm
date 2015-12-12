import Tkinter
import math
import numpy as np
import renderer
import sys

from PIL import Image, ImageTk
from parser import FileParser, NetParser
from renderer import Renderer

class Window(object):
    def __init__(self, fps, renderer):
        self.delay = int(math.ceil(1000. / fps))
        self.renderer = renderer

        p = self.renderer.pixels
        self.root = Tkinter.Tk()
        self.frame = Tkinter.Frame(self.root, width=p, height=p)
        self.frame.pack()
        self.canvas = Tkinter.Canvas(self.frame, width=p, height=p)
        self.canvas.place(x=-2,y=-2)
        self.root.after(0, self.draw)
        self.im = None
        self.photo = None

    def draw(self):
        self.im = Image.fromarray(np.uint8(self.renderer.next()))
        self.photo = ImageTk.PhotoImage(image=self.im)
        self.canvas.create_image(0, 0, image=self.photo, anchor=Tkinter.NW)
        self.root.update()
        self.root.after(self.delay, self.draw)

    def mainloop(self):
        self.root.mainloop()

def main(args):
    if "filename" in args:
        parser = FileParser(args.filename)
    else:
        parser = NetParser(args.host, args.port)
    renderer = Renderer(args.pixels, args.mass, parser)
    window = Window(args.fps, renderer)
    window.mainloop()

if __name__ == "__main__":
    import argparser
    main(argparser.parser().parse_args())
