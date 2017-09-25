---
title: "Problem Description for Diffusion Through Nanopore"
author: Tom Pace
date: September 2017
geometry: left=0.75in, right=0.5in, top=0.5in, bottom=0.5in
---
<!---
This markdown file is intended to be converted to pdf through pandoc with
pandoc --number-sections -o description.pdf description.md
-->

# Background

We seek the effective diffusion constant for a nanoporous membrane.
The pore geometry is variable, as is the diffusion equation.

# Given

## Geometry
We study both a body-centered rectangular lattice of pores,
as well as a face-centered lattice of pores.
The unit cell geometry has two planes of symmetry.
This symmetry is used to reduce the model to only one quarter of the unit cell.

The geometry and geometric variables are shown below.

![Top view of body-centered geometry](./fig_pdf/body-top.pdf){ height=50% }

![Side view of body-centered geometry](./fig_pdf/body-side.pdf){ height=50% }

$\begin{array}{rcl}
Sx & = & 2 Lx =\text{Unit cell x-dimension} \\
Sy & = & 2 Ly =\text{Unit cell y-dimension} \\
Lx & = & \frac{Sx}{2} =\text{Model x-dimension} \\
Ly & = & \frac{Sy}{2} =\text{Model y-dimension} \\
R & = & \text{Pore radius} \\
tm & = & \text{Membrane thickness} = \text{Pore length} \\
H & = & \text{Distance from membrane surface to model boundary}
\end{array}$



_**[[TODO: body-centered side view is incomplete]]**_

_**[[TODO: consider combining views into one image horizontally]]**_

_**[[TODO: face-centered geometry, both views]]**_

## Diffusion Equation

### Unhomogenized Standard Diffusion Equation

### Homogenized Standard Diffusion Equation

# Objective

The effective diffusion constant is found through ...

_**[[TODO]]**_
