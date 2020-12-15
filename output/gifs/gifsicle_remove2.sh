#!/bin/bash
# This script will take an animated GIF and delete every other frame
# Accepts two parameters: input file and output file
# Usage: ./<scriptfilename> input.gif output.gif

# Make a copy of the file
cp $1 $2
echo "copy completed"
# Get the number of frames
numframes=`./gifsicle.exe $1 -I | grep -P "\d+ images" --only-matching | grep -P "\d+" --only-matching`
echo "numframes completed"
# Deletion
./gifsicle.exe "$1" --colors 250 --unoptimize $(seq -f "#%g" 0 2 $numframes) -O2 -o "$2"
#./gifsicle.exe "$1" $(seq -f "#%g" 0 2 $numframes) -o "$2"
echo "script completed"
