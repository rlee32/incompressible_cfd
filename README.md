This repository contains an implementation of J.B. Perot's Exact Fractional 
Step Method, which in contrast to projection methods exactly solves for 
incompressible flows.

This is an old repository of mine in need of revamping. It is written in C.

The code is divided into several parts:

1. STEFS2D: main solver.
2. elgen: grid generator that takes in a geometry and produces a structured 
  grid.
3. con: conversion utility from PLOT3D to STEFS2D format.