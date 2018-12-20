#Introduction
This is a module that allows for generating crack-free isosurface for the block-structured AMR data. We published a paper in IEEE VIS18
For details, plase refer to 
[CPU Isosurface Ray Tracing of Adaptive Mesh Refinement Data](https://ethan0911.github.io/publications/amr)

# How to run the code

Landing Gear Dataset

```bash
#!/bin/bash
IMPI_AMR_DATA=landingGear \
IMPI_AMR_METHOD=octant \
./ospImplicitIsoSurfaceViewer \
/home/sci/feng/Desktop/ws/data/NASA_LandingGear_AMR/cb.osp \
/home/sci/feng/Desktop/ws/data/NASA_LandingGear_AMR/landingGear.obj \
-vp 16.286070 16.446814 0.245150 \
-vu -0.000000 -0.000000 -1.000000 \
-vi 16.430407 16.157639 0.353916 \
"$@"
```

Cosmos Dataset

```bash
#!/bin/bash
IMPI_AMR_DATA=cosmos \
IMPI_AMR_METHOD=octant \
./ospImplicitIsoSurfaceViewer \
/home/sci/feng/Desktop/ws/data/COSMOS_ospray/cosmos.osp \
-vp 260.530029 288.618378 255.794846 \
-vu 1.000000 0.000000 0.000000 \
-vi 248.845215 267.241486 256.334686 \
"$@"
```

Test Dataset

```bash
#!/bin/bash
IMPI_AMR_DATA=chombo \
IMPI_AMR_METHOD=octant \
./ospImplicitIsoSurfaceViewer \
/home/sci/feng/Desktop/ws/data/chombo/chombo_amr.osp \
"$@"
```

`./script.sh` or `mpirun -np <N> ./script.sh --osp:mpi`
