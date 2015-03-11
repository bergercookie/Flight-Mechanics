#!/bin/bash

# Don't know how to write proper Makefiles... 
# So I am writing these little shell scripts t automate the process

echo "************************************************************************"
echo "**Running the Weasel words script**"
bin/weasel *.tex
echo "**Running the Passive voice script**"
bin/passive *.tex
echo "**Running the Duplicates script**"
bin/dups *.tex
echo "************************************************************************"
