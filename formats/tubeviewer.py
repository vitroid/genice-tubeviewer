#!/usr/bin/env python
"""
Develop the sidewalls of dtc-like ice.

Very special utility for dtc study.
"""
from math import pi, atan2
import sys
import itertools as it
import io

import numpy as np
import PIL.ImageDraw as ImageDraw
import PIL.Image as Image
import pairlist as pl

from genice2.decorators import timeit, banner
import genice2.formats


class Format(genice2.formats.Format):
    """
Centers-of-mass of water molecules are output in @AR3A format.
No options available.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {1:self.Hook1}


    @timeit
    @banner
    def Hook1(self, lattice):
        lattice.logger.info("Hook1: Output sidewalls of water nanotubes in dtc-like ice in PNG format.")
        cellmat = lattice.repcell.mat
        zoom = 100

        axes = np.array([(1/8, 1/4), (5/8, 1/4), (3/8, 3/4), (7/8, 3/4)])
        # max distance from the axis to the water monolayer
        maxd = cellmat[0,0] / 4
        radius = maxd - 0.276/2 # approx. radius of the water monolayer
        plusd = maxd + 0.276    # second shell

        grid = pl.determine_grid(cellmat, 0.3)
        pairs = [(i,j) for i,j in pl.pairs_fine(lattice.reppositions, 0.3, cellmat, grid, distance=False)]
        logger.debug(("PAIRS",len(pairs),cellmat,grid))
        sx, sy = int(cellmat[2,2]*zoom),int(4*7*radius*zoom)
        size = (sx//4*4, sy//4*4)
        image = Image.new("RGB", size, '#fff')
        draw  = ImageDraw.Draw(image, "RGBA")

        for j, axis in enumerate(axes):
            dots = dict()
            dots2 = dict()
            for i, rpos in enumerate(lattice.reppositions):
                rd = np.zeros(3)
                rd[:2] = rpos[:2] - axis
                rd -= np.floor(rd+0.5)
                D = rd @ cellmat
                if D@D < maxd**2:
                    a = atan2(D[1], D[0])+3.5
                    z = (rpos @ cellmat)[2]
                    c = a*radius
                    dot = np.array([z,c+j*7*radius])
                    #atom = rpos@cellmat
                    dots[i] = dot
                elif D@D < plusd**2:
                    a = atan2(D[1], D[0])+3.5
                    z = (rpos @ cellmat)[2]
                    c = a*radius
                    dot = np.array([z,c+j*7*radius])
                    #atom = rpos@cellmat
                    dots2[i] = dot
            alldots = {**dots, **dots2}
            for p in dots2:
                r = 0.04
                tl = dots2[p] - r
                br = dots2[p] + r
                tl *= zoom
                br *= zoom
                draw.ellipse([int(x) for x in [tl[0], tl[1], br[0], br[1]]], fill='#f00')
            for p,q in pairs:
                if p in alldots and q in alldots and (p in dots2 or q in dots2):
                    pd = alldots[p] * zoom
                    qd = alldots[q] * zoom
                    D = pd-qd
                    if D@D < zoom**2:
                        draw.line([int(x) for x in [pd[0], pd[1], qd[0], qd[1]]], fill='#f00', width=2)
            for p in dots:
                r = 0.05
                tl = dots[p] - r
                br = dots[p] + r
                tl *= zoom
                br *= zoom
                draw.ellipse([int(x) for x in [tl[0], tl[1], br[0], br[1]]], fill='#000')
            for p,q in pairs:
                if p in dots and q in dots:
                    pd = dots[p] * zoom
                    qd = dots[q] * zoom
                    D = pd-qd
                    if D@D < zoom**2:
                        draw.line([int(x) for x in [pd[0], pd[1], qd[0], qd[1]]], fill='#000', width=3)
        imgByteArr = io.BytesIO()
        tn_image = image.resize((image.width//2, image.height//2), Image.LANCZOS)
        tn_image.save(imgByteArr, format='PNG')
        imgByteArr = imgByteArr.getvalue()
        sys.stdout.buffer.write(imgByteArr)
        lattice.logger.info("Hook1: end.")
