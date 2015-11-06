#!/bin/bash

#Generate geometry without wires
./generate_gdml.pl -w 0 -i microboone-gdml-parameters.xml -o microboone-gdml-fragments.xml
./make_gdml.pl -i microboone-gdml-fragments.xml -o microboonevX_nowires.gdml

#Generate geometry with no wires
./generate_gdml.pl -w 1 -i microboone-gdml-parameters.xml -o microboone-gdml-fragments.xml
./make_gdml.pl -i microboone-gdml-fragments.xml -o microboonevX.gdml

#Copy both back up one level to be seen by LArSoft jobs
cp microboonevX.gdml ..
cp microboonevX_nowires.gdml ..
