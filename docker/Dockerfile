# Runtime for the published GLM-AED executable (glm+ 3.3.1a2, linux/amd64).
# The executable links the netCDF, GD and X11 libraries of Ubuntu 20.04.
# The Intel Fortran runtime libraries are downloaded separately (scripts/01_download.py).
# gcc and libnetcdf-dev compile the state exchange library in native/.
FROM ubuntu:20.04@sha256:8feb4d8ca5354def3d8fce243717141ce31e2c428701f6682bd2fafe15388214
RUN apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      gcc libc6-dev libnetcdf-dev libgd3 libx11-6 \
 && rm -rf /var/lib/apt/lists/*
