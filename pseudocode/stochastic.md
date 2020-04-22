---
title: "Outline of Computational Process for Perturbative Stochastic Homogenization"
# author: Tom Pace
# date: April, 2020
geometry: left=0.75in, right=0.75in, top=0.75in, bottom=0.75in
---
<!---
This markdown file is intended to be converted to pdf through pandoc with
pandoc -o <filename>.pdf <filename>.md
-->

Note: equation references are to Anantharaman and Le Bris, 2011.
However, the notation here is altered.
The primary alterations are:

- Here, we use $D$ instead of $A$ for the diffusion coefficient.
- Here, we use index notation to keep track of vector and matrix components.
- Here, we ascribe variations in $D$ to the problem domain, rather than using different expressions for $D$ to represent different model conditions.

There is also a fair amount of "reading between the lines" required to support this understanding of the process.

The following process is repeated for different supercell sizes, $N$,
with the limit of large $N$ theoretically approaching the correct solution.


# FEM solutions required

1) Obtain the zeroth-order corrector, $w^0_i$.

Solve the usual periodic corrector problem on a single unit cell $X$, with periodic boundary conditions (un-numbered equation between 3.11 and 3.12):

$$-\frac{\partial}{\partial x_k} D^{\text{small}}_{kj}
  \left( \frac{\partial w^0_i}{\partial x_j}
  + \delta_{ij} \right)
  = 0$$

2) Obtain the first-order corrector for supercell size $N$, $w^{1,N}_i$.

Solve the corrector problem on a supercell $X^1_N$ of size $N$ with a single defect, located at the central unit cell, with periodic boundary conditions on the supercell (equation 3.12):

$$-\frac{\partial}{\partial x_k} D^{\text{small}}_{kj}
  \left( \frac{\partial w^{1,N}_i}{\partial x_j}
  + \delta_{ij} \right)
  = 0$$

3) Obtain the second-order correctors, $w^{2,L,N}_i$.

Solve the corrector problem on a supercell $X^{2,L}_N$ of size $N$ with two defects:

- one located at the central unit cell
- one at EACH of the other unit cells, using location index $L$ for the second

with periodic boundary conditions on the supercell (equation 3.13):

$$-\frac{\partial}{\partial x_k} D^{\text{small}}_{kj}
  \left( \frac{\partial w^{2,L,N}_i}{\partial x_j}
  + \delta_{ij} \right)
  = 0$$

(continued...)
\newpage

# Integrated coefficients

1) Obtain the zeroth-order term, $D_{ki}^0$, (equation 3.16 and 2.4).

$$D_{ki}^{0} =
  \frac{1}{|X|}
  \int_X D^{\text{small}}_{kj}
  \left(
    \frac{\partial w^0_i}{\partial x_j} + \delta_{ij} 
  \right) d^d x$$

2) Obtain the first-order term, $D_{ki}^{1,N}$, (equation 3.17).

$$\widetilde{D}_{ki}^{1,N} =
  \frac{1}{|X_N|}
  \int_{X^1_N} D^{\text{small}}_{kj}
  \left(
    \frac{\partial w^{1,N}_i}{\partial x_j} + \delta_{ij}
  \right) d^d x$$

$$D_{ki}^{1,N} = \widetilde{D}_{ki}^{1,N} - D_{ki}^0$$

3) Obtain the second-order term, $D_{ki}^{2,N}$, (equation 3.18).

$$\widetilde{D}_{ki}^{2,L,N} =
  \frac{1}{|X_N|}
  \int_{X^{2,L}_N} D^{\text{small}}_{kj}
  \left(
    \frac{\partial w^{2,L,N}_i}{\partial x_j} + \delta_{ij}
  \right) d^d x$$

$$\overline{\widetilde{D}}_{ki}^{2,N} = 
\frac{1}{N^d - 1} \sum_{L=1,L\neq\text{center}}^{N^d} \left( \widetilde{D}_{ki}^{2,L,N} \right)$$

$$D_{ki}^{2,N} = \frac{1}{2} \left(
  \overline{\widetilde{D}}_{ki}^{2,N}
  - 2 \widetilde{D}_{ki}^{1,N} + D_{ki}^0
  \right)$$


# Final result

Combine the terms of the integrated coefficients (equation 3.15).

$$D_{ki}^{\text{large}} = \lim_{N\to\infty}
\left(
  D_{ki}^0 + \eta D_{ki}^{1,N} + \eta^2 D_{ki}^{2,N}
\right)$$

Where $\eta$ is the defect probability.

# Notes

When expressed in these forms,
it is clear that the form of the corrector problem is the same at each order,
and the same integral is also required for each order
(the $\widetilde{D}$ values).

Accordingly, the following FEM simulations are performed:

- One simulation of a non-defect single unit cell
- One simulation of a single-defect supercell, with the defect at the center
- $N^d-1$ simulations of two-defect supercells, with one defect at the center and the other elsewhere

For each of these simulations,
the corrector problem is the same as the one for periodic homogenization,
and the usual homogenization integral is evaluated as well.
This yields all the necessary $\widetilde{D}$ values.

The final result can then be computed from these simulation results
with simple arithmetic:

- The average of the two-defect results is computed, $\overline{\widetilde{D}}_{ki}^{2,N}$.
- $D_{ki}^{1,N}$ and $D_{ki}^{2,N}$ are computed.
- Finally, $D_{ki}^{\text{large}}$ is computed.
