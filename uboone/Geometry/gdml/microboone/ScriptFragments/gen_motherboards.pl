#!/usr/bin/perl

print MOTHERBOARDS <<EOF;

<?xml version='1.0'?>
<gdml>
<solids>
 <box name="ElectroBoard"
  lunit="cm"
  z="30"
  x="18"
  y=".5"/>
 <box name="ElectroSetPiece"
  lunit="cm"
  z="1"
  x="15"
  y="5"/>
 <box name="ElectroSet"
  lunit="cm"
  z="30"
  x="15"
  y="5"/>
 <box name="Nub1"
  lunit="cm"
  x="1"
  y="1.5"
  z=".5"/>
 <box name="Nub2"
  lunit="cm"
  x="2"
  y="1.5"
  z=".5"/>
 <box name="Nub3"
  lunit="cm"
  x="3"
  y="1.5"
  z=".5"/>
 <box name="Nub4"
  lunit="cm"
  x="1.5"
  y="1.5"
  z=".5"/>
 <box name="ElectroCard"
  lunit="cm"
  x="9"
  y="3"
  z=".5"/>
 <box name="Electronics"
  lunit="cm"
  x="15"
  y="10"
  z="800"/>
<xtru name="ElectroBar" lunit="cm" >
  <twoDimVertex x="1.3" y="1.9" />
  <twoDimVertex x="1.3" y="2.1" />
  <twoDimVertex x="-1.3" y="2.1" />
  <twoDimVertex x="-1.3" y="-2.1" />
  <twoDimVertex x="1.3" y="-2.1" />
  <twoDimVertex x="1.3" y="-1.9" />
  <twoDimVertex x="-1.1" y="-1.9" />
  <twoDimVertex x="-1.1" y="1.9" />
  <section zOrder="0" zPosition="-400" xOffset="0" yOffset="0" scalingFactor="1" />
  <section zOrder="1" zPosition="400" xOffset="0" yOffset="0" scalingFactor="1" />
</xtru>
</solids>
<structure>
 <volume name="volElectroBar">
  <materialref ref="G10"/>
  <solidref ref="ElectroBar"/>
 </volume>
 <volume name="volElectroBoard">
  <materialref ref="SILICON_Si"/>
  <solidref ref="ElectroBoard"/>
 </volume>
 <volume name="volNub1">
  <materialref ref="G10"/>
  <solidref ref="Nub1"/>
 </volume>
 <volume name="volNub2">
  <materialref ref="G10"/>
  <solidref ref="Nub2"/>
 </volume>
 <volume name="volNub3">
  <materialref ref="G10"/>
  <solidref ref="Nub3"/>
 </volume>
 <volume name="volNub4">
  <materialref ref="G10"/>
  <solidref ref="Nub4"/>
 </volume>
 <volume name="volElectroCard">
  <materialref ref="SILICON_Si"/>
  <solidref ref="ElectroCard"/>
 </volume>
 <volume name="volElectroSetPiece">
  <materialref ref="LAr"/>
  <solidref ref="ElectroSetPiece"/>
  <physvol>
   <volumeref ref="volNub1"/>
   <position name="posNub1" unit="cm" x="6.5" y="-.75" z="-.25"/>
  </physvol>
  <physvol>
   <volumeref ref="volNub1"/>
   <position name="posNub1" unit="cm" x="4.5" y="-.75" z="-.25"/>
  </physvol>
  <physvol>
   <volumeref ref="volNub2"/>
   <position name="posNub2" unit="cm" x="1.5" y="-.75" z="-.25"/>
  </physvol>
  <physvol>
   <volumeref ref="volNub3"/>
   <position name="posNub3" unit="cm" x="-1.5" y="-.75" z="-.25"/>
  </physvol>
  <physvol>
   <volumeref ref="volNub4"/>
   <position name="posNub4" unit="cm" x="-4.5" y="-.75" z="-.25"/>
  </physvol>
  <physvol>
   <volumeref ref="volElectroCard"/>
   <position name="posElectroCard" unit="cm" x="-1.5" y="1" z=".25"/>
  </physvol>
 </volume> 
 <volume name="volElectroSet">
  <materialref ref="LAr"/>
  <solidref ref="ElectroSet"/>
  <physvol>
   <volumeref ref="volElectroBoard"/>
   <position name="posElectroBoard" unit="cm" x="0" y="-1.5" z="0"/>
  </physvol>
EOF

$count=0;
while ( $count < 10 ) {
    
    $margin=$count*1.5;

    print MOTHERBOARDS <<EOF;
  <physvol> 
    <volumeref ref="volElectroSetPiece"/> 
    <position name="posElectroSetPieceP$count" unit="cm" x="0" y="0" z="$margin" /> 
  </physvol> 
EOF

  if ( $count != 0 ) {
        print MOTHERBOARDS <<EOF;
   <physvol> 
     <volumeref ref="volElectroSetPiece"/> 
     <position name="posElectroSetPieceN$count" unit="cm" x="0" y="0" z="-$margin" /> 
   </physvol> 
EOF
  }
  $count++;
}

    print MOTHERBOARDS <<EOF
 </volume>
 <volume name="volElectronics">
  <materialref ref="LAr"/>
  <solidref ref="Electronics"/>
  <physvol>
   <volumeref ref="volElectroBar"/>
   <position name="posElectroBar1" unit="cm" x="6.5" y="-2" z="0"/>
  </physvol>
 </volume> 
</structure>
</gdml>
EOF

$count=0;
while ( $count < 9 ) {
  $margin=$count*32;
    print MOTHERBOARDS <<EOF;
  <physvol> 
   <volumeref ref="volElectroSet"/> 
   <position name="posElectroSetP$count" unit="cm" x="0" y="2" z="$margin"/> 
  </physvol> 
EOF
  if ( $count != 0 ) {
        print MOTHERBOARDS <<EOF;
    <physvol> 
     <volumeref ref="volElectroSet"/> 
     <position name="posElectroSetN$count" unit="cm" x="0" y="2" z=."$margin"/> 
    </physvol> 
EOF
  }
  $count++;
}
