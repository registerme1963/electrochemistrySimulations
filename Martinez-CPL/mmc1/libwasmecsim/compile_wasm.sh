#!/bin/bash

em++ main.cpp system.cpp electrodes.cpp experiment.cpp coefs_alpha_beta.cpp simulation.cpp -o libwasmecsim.js -O3 -std=c++11 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['ccall', 'cwrap']"
