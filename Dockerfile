# Dockerfile — canonical analysis environment (reorg plan v2, phase-6 unit)
#
# Pinning (single source: reorg plan v2 appendix C; D-001/002/003/018/019 are
# the historical basis):
#   base   : ubuntu:noble-20260410 (date tag = immutable build; digest at
#            adoption, 2026-08-09:
#            sha256:c4a8d5503dfb2a3eb8ab5f807da5bc69a85730fb49b5cfca2330194ebcc41c7b)
#   apt    : snapshot.ubuntu.com/20260410T000000Z, Level 2 (new installs only,
#            no upgrade), matching the base image's apt-layer build date
#   R      : 4.5.3 built from the CRAN source tarball (sha256-verified below),
#            --without-recommended-packages; recommended set comes from P3M
#            like everything else, so one mechanism pins every package
#   R pkgs : P3M date snapshot 2026-04-09 (CRAN + Bioconductor 3.22),
#            baked into the image (no host-mounted library)
#   BLAS   : reference BLAS/LAPACK (libblas-dev/liblapack-dev from the apt
#            snapshot). OpenBLAS is ABSENT BY DESIGN: amd64 bit-identity is
#            then a property of the build, not of runtime dispatch
#            (qualification measurements: reorg plan v2 B.12 — zero time cost,
#            all persisted artifacts bit-identical, MUREN BLAS-free).
#            The build fails if any OpenBLAS package is present or if R links
#            it.
#
# Usage: the repository is mounted at runtime, never copied into the image:
#   docker build -t rebc-r453:refblas .
#   docker run -d --name rebc-r453 -v /path/to/rstudio-docker:/workspace \
#     rebc-r453:refblas sleep infinity
#   docker exec rebc-r453 Rscript scripts/<step>.R

FROM ubuntu:noble-20260410

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 DEBIAN_FRONTEND=noninteractive

# --- apt: pin to the snapshot matching the base image build date (Level 2) ---
RUN printf '%s\n' \
      'Types: deb' \
      'URIs: https://snapshot.ubuntu.com/ubuntu/20260410T000000Z' \
      'Suites: noble noble-updates noble-backports noble-security' \
      'Components: main restricted universe multiverse' \
      'Signed-By: /usr/share/keyrings/ubuntu-archive-keyring.gpg' \
      > /etc/apt/sources.list.d/ubuntu.sources

# Bootstrap: the snapshot host is HTTPS-only and the base image carries no CA
# store, so ca-certificates itself is fetched once with peer verification off
# (recorded deviation, same as the 4-1b provisioning; every later fetch is
# verified).
RUN apt-get update -o Acquire::https::Verify-Peer=false \
 && apt-get install -o Acquire::https::Verify-Peer=false \
      -y --no-install-recommends ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# --- toolchain, R build dependencies, and the reference BLAS/LAPACK ---------
# Same set the dev container (4-1b) proved out, minus every OpenBLAS package,
# plus libblas-dev/liblapack-dev (reference implementations, 3.12.0 at this
# snapshot — the exact libraries the B.12 qualification validated).
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential gfortran wget xz-utils cmake \
      libreadline-dev libpcre2-dev zlib1g-dev libbz2-dev liblzma-dev \
      libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev libncurses-dev \
      libcairo2-dev libfontconfig1-dev libfreetype6-dev \
      libpng-dev libjpeg-dev libtiff-dev \
      libblas-dev liblapack-dev \
 && rm -rf /var/lib/apt/lists/* \
 && if dpkg -l | grep -qi openblas; then \
      echo 'OpenBLAS must not be present in this image' >&2; exit 1; fi

# --- R 4.5.3 from CRAN source, linked against the reference BLAS/LAPACK ------
ARG R_VERSION=4.5.3
ARG R_SHA256=aa5c1ed4293c7271ac513d654670356ac0e8a6ad5e42be014365d11150b5b8f2
RUN cd /tmp \
 && wget -q "https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz" \
 && echo "${R_SHA256}  R-${R_VERSION}.tar.gz" | sha256sum -c - \
 && tar xzf "R-${R_VERSION}.tar.gz" \
 && cd "R-${R_VERSION}" \
 && ./configure --prefix=/usr/local \
      --enable-R-shlib --with-blas --with-lapack \
      --without-recommended-packages --without-x \
 && make -j"$(nproc)" \
 && make install \
 && cd / && rm -rf /tmp/R-${R_VERSION} /tmp/R-${R_VERSION}.tar.gz

# --- package pinning: P3M date snapshot 2026-04-09 (CRAN + Bioc 3.22) --------
RUN printf '%s\n' \
      '# P3M snapshot 2026-04-09 (reorg plan v2 appendix C)' \
      'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2026-04-09"))' \
      'options(BioC_mirror = "https://packagemanager.posit.co/bioconductor/2026-04-09")' \
      'options(BIOCONDUCTOR_CONFIG_FILE = "https://packagemanager.posit.co/bioconductor/2026-04-09/config.yaml")' \
      'Sys.setenv("R_BIOC_VERSION" = "3.22")' \
      > /usr/local/lib/R/etc/Rprofile.site

# --- R packages: explicit set from the active pipeline code, baked in --------
COPY docker/install_packages.R /tmp/docker-build/
RUN mkdir -p /usr/local/lib/R/site-library \
 && Rscript /tmp/docker-build/install_packages.R

# --- build-time verification: versions, BLAS linkage, loadability, Rcpp ------
# (separate COPY so verification changes do not re-trigger the install layer)
COPY docker/verify_environment.R docker/versions.tsv /tmp/docker-build/
RUN Rscript /tmp/docker-build/verify_environment.R \
 && rm -rf /tmp/docker-build

# --- PID-1 reaper + figure font ---------------------------------------------
# mclapply forks orphan at script exit; without a reaping init they accumulate
# as zombies in long-lived containers (917 observed in the dev container after
# ~1 day). tini as ENTRYPOINT reaps them regardless of how the container is
# invoked. (Late layer on purpose: keeps the R build cache during iteration.)
# fonts-liberation: Arial-metric sans-serif for figure text (Springer Nature
# artwork guide: Helvetica/Arial, 5-7 pt). Display layer only; same pinned
# apt snapshot (added 2026-08-28, researcher Go).
RUN apt-get update \
 && apt-get install -y --no-install-recommends tini fonts-liberation \
 && rm -rf /var/lib/apt/lists/*

USER ubuntu
WORKDIR /workspace
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]
