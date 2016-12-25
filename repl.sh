#!/bin/bash
#SBATCH -J repl
#SBATCH -o repl.Rout
#SBATCH -e repl.Rout
#SBATCH -p parallel-short

R --save  <./nobuild/jss_revision/replication.R
