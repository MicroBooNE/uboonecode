#!/usr/bin/perl

gen_wheelsupport();

sub gen_wheelsupport 
{

$WHEELS = "micro-wheelsupport.gdml" . $suffix . ".gdml";
push (@gdmlFiles, $WHEELS); # Add file to list of WHEELS fragments
$WHEELS = ">" . $WHEELS;
open(WHEELS) or die("Could not open file $WHEELS for writing");

print WHEELS <<EOF;


<?xml version='1.0'?>
<gdml>
<solids>
 <box name="WheelSupport"
  x="240"
  y="30"
  z="16"
  lunit="cm"/>
 <xtru name="WheelSupportPiece1" lunit="cm" >	
  <twoDimVertex x="89.15" y="0" />
  <twoDimVertex x="119.7" y="10.9" />
  <twoDimVertex x="-114.8" y="10.9" />
  <twoDimVertex x="-84.25" y="0" />
  <section zOrder="0" zPosition="-7.8/2" xOffset="0" yOffset="0" scalingFactor="1" />
  <section zOrder="1" zPosition="7.8/2" xOffset="0" yOffset="0" scalingFactor="1" />
 </xtru>
 <box name="WheelSupportPiece2" 
  x="32"
  y="5"
  z="1"
  lunit="cm"/>
 <box name="WheelSupportPiece3" 
  x="16"
  y="1"
  z="7"
  lunit="cm"/>
 <xtru name="WheelSupportPiece4" lunit="cm" >	
  <twoDimVertex x="3.5" y="0" />
  <twoDimVertex x="3.5" y="10" />
  <twoDimVertex x="8" y="10" />
  <twoDimVertex x="8" y="-14.8" />
  <twoDimVertex x="-8" y="-14.8" />
  <twoDimVertex x="-8" y="10" />
  <twoDimVertex x="-3.5" y="10" />
  <twoDimVertex x="-3.5" y="0" />
  <section zOrder="0" zPosition="-7.8/2" xOffset="0" yOffset="0" scalingFactor="1" />
  <section zOrder="1" zPosition="7.8/2" xOffset="0" yOffset="0" scalingFactor="1" />
 </xtru>
</solids>

<structure>
 <volume name="volWheelSupportPiece1">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="WheelSupportPiece1"/>
 </volume>
 <volume name="volWheelSupportPiece2">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="WheelSupportPiece2"/>
 </volume>
 <volume name="volWheelSupportPiece3">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="WheelSupportPiece3"/>
 </volume>
 <volume name="volWheelSupportPiece4">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="WheelSupportPiece4"/>
 </volume>
 <volume name="volWheelSupport">
  <materialref ref="LAr"/>
  <solidref ref="WheelSupport"/>
  <physvol>
   <volumeref ref="volWheelSupport1"/>
   <position name="posWheelSupport1" unit="cm" x="0" y="0" z="0"/>
   <rotationref ref="rPlus90AboutY" />
  </physvol>
  <physvol>
   <volumeref ref="volWheelSupport2"/>
   <position name="posWheelSupport2a" unit="cm" x="12+79.7" y="13.4" z="3"/>
  </physvol>
  <physvol>
   <volumeref ref="volWheelSupport2"/>
   <position name="posWheelSupport2b" unit="cm" x="12+79.7" y="13.4" z="-3"/>
  </physvol>
  <physvol>
   <volumeref ref="volWheelSupport2"/>
   <position name="posWheelSupport2c" unit="cm" x="-(12+74.8)" y="13.4" z="3"/>
  </physvol>
  <physvol>
   <volumeref ref="volWheelSupport2"/>
   <position name="posWheelSupport2d" unit="cm" x="-(12+74.8)" y="13.4" z="-3"/>
  </physvol>
  <physvol>
   <volumeref ref="volWheelSupport3"/>
   <position name="posWheelSupport3a" unit="cm" x="(8+79.7)" y="16.4" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volWheelSupport3"/>
   <position name="posWheelSupport3b" unit="cm" x="-(8+74.8)" y="16.4" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volWheelSupport4"/>
   <position name="posWheelSupport4" unit="cm" x="-52.15" y="0" z="0"/>
   <rotationref ref="rPlus90AboutY" />
  </physvol>
  <physvol>
   <volumeref ref="volWheelSupport4"/>
   <position name="posWheelSupport4" unit="cm" x="-52.15" y="0" z="0"/>
   <rotationref ref="rPlus90AboutY" />
  </physvol>
 </volume>
</structure>
</gdml>
EOF

close(WHEELS);
}
