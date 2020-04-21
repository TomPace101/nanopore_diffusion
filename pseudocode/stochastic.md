---
title: "Pseudocode for Stochastic Homogenization"
# author: Tom Pace
# date: April, 2020
geometry: left=0.75in, right=0.5in, top=0.5in, bottom=0.5in
---
<!---
This markdown file is intended to be converted to pdf through pandoc with
pandoc -o <filename>.pdf <filename>.md
-->

Note: equation references are to Anantharaman and Le Bris, 2011.

# FEM solutions required

1) Obtain the zeroth-order corrector, $w^0_i$.

Solve the usual periodic corrector problem on a single unit cell, with periodic boundary conditions (un-numbered equation between 3.11 and 3.12):

$$-\frac{\partial}{\partial x_k} D^{\text{small}}_{kj}
  \left( \frac{\partial w^0_i}{\partial x_j}
  + \delta_{ij} \right)
  = 0$$

2) Obtain the first-order corrector for supercell size $N$, $w^{1,N}_i$.

Solve the corrector problem on a supercell of size $N$ with a single defect, located at the central unit cell, with periodic boundary conditions on the supercell (equation 3.12):

$$-\frac{\partial}{\partial x_k} D^{\text{small}}_{kj}
  \left( \frac{\partial w^{1,N}_i}{\partial x_j}
  + \delta_{ij} \right)
  = 0$$

3) Obtain the second-order correctors for supercell size $N$, $w^{2,L,N}_i$.

Solve the corrector problem on a supercell of size $N$ with two defects:

- one located at the central unit cell
- one at EACH of the other unit cells, using location index $L$ for the second

with periodic boundary conditions on the supercell (equation 3.13):

$$-\frac{\partial}{\partial x_k} D^{\text{small}}_{kj}
  \left( \frac{\partial w^{2,L,N}_i}{\partial x_j}
  + \delta_{ij} \right)
  = 0$$

# Integrated coefficients

