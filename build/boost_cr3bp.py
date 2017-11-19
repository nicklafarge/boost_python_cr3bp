#!/usr/bin/env python

import boost_cr3bp

example_ic = [-0.27, -0.42, 0, 0.3, -1.0, 0]
example_t = [0, 1]
mu = 0.0121505842699

cr3bp = boost_cr3bp.CR3BP()
example_traj = cr3bp.prop(example_ic, example_t, mu, 6, 2, 1e-12, 1e-5)
print example_traj