#!/bin/bash

cmake -DGeant4_DIR= ~/Downloads/Geant4-10.1.0-Linux/lib64 $1

make -j4

# useless comment
