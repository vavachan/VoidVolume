# VoidVolume
C++ code to compute the volume of an overlapping set of sphere. Equivalently, it can compute the void volume. 

Ref : **Sastry, Srikanth, et al. "Statistical geometry of particle packings. I. Algorithm for exact determination of connectivity, volume, and surface areas of void space in monodisperse and polydisperse sphere packings." Physical Review E 56.5 (1997): 5524.**

```
g++ voidVolume3D.cpp -std=c++11 -O3

./a.out [filename] [probe_radius]
```

Code can handle periodic boundary conditions, but one has to ensure that any sphere is not forming a Delanuay tessellation with another particle and its image as well. 

Use VMD to visualize the Delanuay Tessellations. 
