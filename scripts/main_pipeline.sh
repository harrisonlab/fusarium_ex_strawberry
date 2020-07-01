#!/bin/bash

export PATH=${PATH}:~/git_repos/assembly_fusarium_ex/scripts

# spades assembly step
spades_step.sh

# quast assembly QC step
quast_step.sh
