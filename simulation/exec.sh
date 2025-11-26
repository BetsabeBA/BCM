#!/bin/bash
# Example batch script to run BCM_simulation.R and log output

Rscript simulation/BCM_simulation.R > simulation/log/sim_log.txt 2>&1
