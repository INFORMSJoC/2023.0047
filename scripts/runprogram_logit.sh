#!/bin/sh
homepath="/home/user"

# Path to Julia binary and to script file
executable="$homepath/programs/julia/julia-1.1.0/bin/julia"
scriptfile="$homepath/path/to/main_file/mainLogit.jl"

# set ld path
export LD_LIBRARY_PATH="$homepath/programs/julia/julia-1.1.0/lib":$LD_LIBRARY_PATH

# get parameters from logit_run.sh
params=$1


# execute program
$executable $scriptfile $params



