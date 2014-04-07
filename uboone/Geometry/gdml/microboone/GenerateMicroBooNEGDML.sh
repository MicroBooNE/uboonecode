#!/bin/bash

#Generate geometry without wires
./generate_gdml.pl -w 0 -i microboone-gdml-parameters.xml -o microboone-gdml-fragments.xml
./make_gdml.pl -i microboone-gdml-fragments.xml -o microboone_nowires.gdml

#Generate geometry with no wires
./generate_gdml.pl -w 1 -i microboone-gdml-parameters.xml -o microboone-gdml-fragments.xml
./make_gdml.pl -i microboone-gdml-fragments.xml -o microboone.gdml

#Generate geometry with no cryostat 
./generate_gdml.pl -w 0 -c 0 -i microboone-gdml-parameters.xml -o microboone-gdml-fragments.xml
./make_gdml.pl -i microboone-gdml-fragments.xml -o microboone_LArTF.gdml

#Copy both back up one level to be seen by LArSoft jobs
cp microboone.gdml ..
cp microboone_nowires.gdml ..
cp microboone_LArTF.gdml ..
