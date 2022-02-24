#!/bin/bash
rocprof -o results.${SLURM_PROCID}.csv --obj-tracking on --hip-trace --parallel-kernels $@
