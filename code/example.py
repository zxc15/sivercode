:#!/usr/bin/env python

import lhapdf

p = lhapdf.mkPDF("CT14lo",0)
print(p.alphasQ(0.1))
print(p.xfxQ(21, 1e-3, 1e4))
print(p.flavors())


