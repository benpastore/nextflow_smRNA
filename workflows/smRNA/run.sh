#!/bin/bash

set -e

function usage
{
    echo "usage: sh run.sh -d -s SOME_MORE_ARGS [-y YET_MORE_ARGS || -h]"
    echo "   ";
    echo "  -h | --help : This message";
}

function parse_args
{
  # positional args
  args=()

  # named args
  while [ "$1" != "" ]; do
    case "$1" in
      -h | --help )         usage;            exit;; # quit and show usage
      * )                   args+=("$1")             # if no match, add it to the positional args
    esac
    shift 
    # move to next kv pair
  done

  # restore positional args
  set -- "${args[@]}"
}

function run
{

parse_args "$@"
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo """#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow $SCRIPT_DIR/main.nf $@


""" > $SCRIPT_DIR/submit.sbatch
sbatch $SCRIPT_DIR/submit.sbatch

}
  
run "$@";