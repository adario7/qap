#!/bin/bash
#SBATCH --job-name QAP
#SBATCH --partition arrow 
#SBATCH --ntasks 1
#SBATCH --mem 14GB
#SBATCH --time 12:05:00

set -x

# warm up processors
sudo cpupower frequency-set -g performance

sleep 0.1

stress-ng -c 4 --cpu-ops=100

# set limits
ulimit -v 14884864



# Generate name by joining all arguments with an underscore
name=$(echo "$@" | tr ' ' '_')
# Extract command and input file
command="$1"
input_file="$2"
# Extract all parameters after the second one
shift 2
params="$@"
# Run the command with the provided parameters
./src/build/$command $params < "instance/$input_file" > "out/o_$name" 2> "out/e_$name"



# back to power saving mode
sudo cpupower frequency-set -g powersave
