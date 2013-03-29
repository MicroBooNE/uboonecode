#!/usr/bin/perl

    $field_cage_width		=	230;
    $field_cage_length		=	1050;
    $field_cage_loop_interval	=	4; 	# =1 is normal, =4 skips 3/4
    $spacers_on_off		= 	"off"; 	# "on" or "off" for tube spacers (off saves time)

    # Set up the output file.
    $FIELDCAGE = "micro-fieldcage.gdml";
    #push (@defFiles, $FIELDCAGE); # Add file to list of constant files
    $FIELDCAGE = ">" . $FIELDCAGE;
    open(FIELDCAGE) or die("Could not open file $FIELDCAGE for writing");

    # Print the Field Cage constants
    print FIELDCAGE <<EOF;
<define>
 <constant name="kFieldCageTPCClearance"    	value="5*kInch" />

 <constant name="kFieldCageTubeRadius"  	value="0.5*kInch" />
 <constant name="kFieldCageTubeThickness" 	value="0.25*kInch" />
 <constant name="kFieldCageBeamDepth"           value="12.5"/>
 <constant name="kFieldCageBeamWidth"           value="2.5"/>
 <constant name="kFieldCageCrossDepth"          value="2.5"/>
 <constant name="kFieldCageCrossWidth"          value="4"/>
 <constant name="kFieldCageCrossLength"         value="308"/>

 <constant name="kTPCTotalLength"         	value="1050"/>
 <constant name="kTPCTotalWidth"         	value="240"/>
 <constant name="kTPCTotalHeight"         	value="220"/>

 <constant name="kFieldCageLoopLength"       	value="kTPCTotalLength+2*(kFieldCageTPCClearance+2*kFieldCageTubeRadius)"/>
 <constant name="kFieldCageLoopWidth"   	value="kTPCTotalWidth"/>
 <constant name="kFieldCageLoopHeight"       	value="kTPCTotalHeight+2*(kFieldCageTPCClearance+2*kFieldCageTubeRadius)"/>

 <constant name="kFieldCageCornerRadius"	value="0.5*kFieldCageTPCClearance"/>
 <constant name="kFieldCageCornerY"       	value="(0.5*kFieldCageLoopHeight)-kFieldCageCornerRadius-kFieldCageTubeRadius"/>
 <constant name="kFieldCageCornerZ"       	value="(0.5*kFieldCageLoopLength)-kFieldCageCornerRadius-kFieldCageTubeRadius"/>

 <constant name="kFieldCageHeight"       	value="kFieldCageLoopHeight+2*(kFieldCageBeamDepth+kFieldCageCrossDepth)"/>
 <constant name="kFieldCageLength"       	value="kFieldCageLoopLength+2*(kFieldCageBeamDepth+kFieldCageCrossDepth)"/>
 <constant name="kFieldCageWidth"      		value="kFieldCageLoopWidth"/>

 <constant name="kFieldCageBeamYInt"            value="0.5*(kFieldCageLoopHeight-50)"/>
 <constant name="kFieldCageBeamZPos"            value="0.5*(kFieldCageLoopLength)"/>
 <constant name="kFieldCageBeamYPos"            value="0.5*(kFieldCageLoopHeight)"/>
 <constant name="kFieldCageBeamZInt"            value="0.5*(kFieldCageLoopLength-kFieldCageCornerRadius-kFieldCageTubeRadius)"/>

 <constant name="kFieldCageCrossYPos"           value="0.5*(kFieldCageLoopHeight+kFieldCageCrossDepth)+kFieldCageBeamDepth"/>
 <constant name="kFieldCageCrossZPos"           value="0.5*(kFieldCageLoopLength+kFieldCageCrossDepth)+kFieldCageBeamDepth"/>
</define>
EOF

    # Prints Field Cage solids
    print FIELDCAGE <<EOF;
<solids>
 <box name="FieldCage"
  x="kFieldCageWidth"
  y="kFieldCageHeight"
  z="kFieldCageLength"
  lunit="cm"/>
 <box name="FieldCageLoop"
  x="2*kFieldCageTubeRadius"
  y="kFieldCageLoopHeight"
  z="kFieldCageLoopLength"
  lunit="cm"/>
 <tube name="FieldCageTubeZ"
  rmin="kFieldCageTubeRadius-kFieldCageTubeThickness"
  rmax="kFieldCageTubeRadius"  
  z="kFieldCageLoopLength-2*(kFieldCageCornerRadius+kFieldCageTubeRadius)" 
  deltaphi="2*kPi" 
  aunit="rad" 
  lunit="cm"/> 
 <tube name="FieldCageTubeY" 
  rmin="kFieldCageTubeRadius-kFieldCageTubeThickness"
  rmax="kFieldCageTubeRadius"  
  z="kFieldCageLoopHeight-2*(kFieldCageCornerRadius+kFieldCageTubeRadius)"  
  deltaphi="2*kPi" 
  aunit="rad" 
  lunit="cm"/> 
 <torus name="FieldCageCorner"
  rmin="kFieldCageTubeRadius-kFieldCageTubeThickness"
  rmax="kFieldCageTubeRadius"  
  rtor="kFieldCageCornerRadius"
  deltaphi="0.5*kPi"
  aunit="rad"
  lunit="cm"/>
 <box name="FieldCageBarFlat"
  x="kFieldCageLoopWidth"
  y="kFieldCageBeamDepth"
  z="kFieldCageBeamWidth"
  lunit="cm"/>
 <box name="FieldCageBarLoopSpacer"
  x="4-2*kFieldCageTubeRadius"
  y="2*kFieldCageTubeRadius"
  z="kFieldCageBeamWidth"
  lunit="cm"/>
 <box name="FieldCageBarLoopSpacers"
  x="kFieldCageWidth"
  y="2*kFieldCageTubeRadius"
  z="kFieldCageBeamWidth"
  lunit="cm"/>
 <box name="FieldCageBarFlat2"
  x="kFieldCageLoopWidth"
  y="0.5*kFieldCageBeamDepth"
  z="kFieldCageBeamWidth"
  lunit="cm"/>
 <box name="FieldCageBarAngle"
  x="kFieldCageCrossLength"
  y="kFieldCageCrossWidth"
  z="kFieldCageCrossDepth"
  lunit="cm"/>
 <box name="FieldCageCrossBox"
  x="36"
  y="30"
  z="kFieldCageCrossDepth"
  lunit="cm"/>
 <box name="FieldCageCross"
  x="kFieldCageLoopWidth"
  y="kFieldCageLoopHeight"
  z="kFieldCageCrossDepth"
  lunit="cm"/>

 <xtru name="Triangle1" lunit="cm" >
  <twoDimVertex x="0" y="0" />
  <twoDimVertex x="-15" y="0" />
  <twoDimVertex x="-15" y="-15" />
  <section zOrder="0" zPosition="0" xOffset="0" yOffset="0" scalingFactor="1" />
  <section zOrder="1" zPosition="2" xOffset="0" yOffset="0" scalingFactor="1" />
 </xtru>
 <xtru name="Triangle2" lunit="cm" >
  <twoDimVertex x="0" y="0" />
  <twoDimVertex x="15" y="0" />
  <twoDimVertex x="15" y="-15" />
  <section zOrder="0" zPosition="0" xOffset="0" yOffset="0" scalingFactor="1" />
  <section zOrder="1" zPosition="2" xOffset="0" yOffset="0" scalingFactor="1" />
 </xtru>
</solids> 
EOF

    # Prints Field Cage tube loop sub-structure
    print FIELDCAGE <<EOF;
<structure> 
 <volume name="volFieldCageTubeTop"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeZ"/> 
 </volume> 
 <volume name="volFieldCageTubeBot"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeZ"/> 
 </volume> 
 <volume name="volFieldCageTubeFront"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeY"/> 
 </volume> 
 <volume name="volFieldCageTubeBack"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageTubeY"/> 
 </volume> 
 <volume name="volFieldCageCorner"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="FieldCageCorner"/> 
 </volume> 

 <volume name="volFieldCageLoop">
  <materialref ref="LAr"/> 
  <solidref ref="FieldCageLoop"/> 
  <physvol>
   <volumeref ref="volFieldCageTubeTop"/>
   <position name="posFieldCageTubeTop" unit="cm" x="0" y="0.5*(kFieldCageLoopHeight-2*kFieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBot"/>
   <position name="posFieldCageTubeBot" unit="cm" x="0" y="-0.5*(kFieldCageLoopHeight-2*kFieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeFront"/>
   <position name="posFieldCageTubeFront" unit="cm" x="0" y="0" z="0.5*(kFieldCageLoopLength-2*kFieldCageTubeRadius)"/>
   <rotation name="rFieldCageVert" aunit="rad" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBack"/>
   <position name="posFieldCageTubeBack" unit="cm" x="0" y="0" z="-0.5*(kFieldCageLoopLength-2*kFieldCageTubeRadius)"/>
   <rotation name="rFieldCageVert" aunit="rad" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCorner"/>
   <position name="posFieldCageCornerBackTop" unit="cm" x="0" y="kFieldCageCornerY" z="-kFieldCageCornerZ"/>
   <rotation name="rFieldCageCorner" aunit="rad" x="90" y="-90" z="-90"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCorner"/>
   <position name="posFieldCageCornerBackBot" unit="cm" x="0" y="-kFieldCageCornerY" z="-kFieldCageCornerZ"/>
   <rotationref ref="rMinus90AboutYPlus90AboutZ"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCorner"/>
   <position name="posFieldCageCornerFrontTop" unit="cm" x="0" y="kFieldCageCornerY" z="kFieldCageCornerZ"/>
   <rotationref ref="rPlus90AboutY"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCorner"/>
   <position name="posFieldCageCornerFrontBot" unit="cm" x="0" y="-kFieldCageCornerY" z="kFieldCageCornerZ"/>
   <rotationref ref="rMinus90AboutYPlus180AboutZ"/>
  </physvol>
 </volume> 
EOF

    # Prints Field Cage bar definitions
    print FIELDCAGE <<EOF;
 <volume name="volFieldCageBarAngle1">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="FieldCageBarAngle"/> 
 </volume>  
 <volume name="volFieldCageBarAngle2">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="FieldCageBarAngle"/>
 </volume>
 <volume name="volFieldCageCrossBox">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="FieldCageCrossBox"/>
 </volume>
 <volume name="volTriangle1"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="Triangle1"/> 
 </volume> 
 <volume name="volTriangle2"> 
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/> 
  <solidref ref="Triangle2"/> 
 </volume> 

 <volume name="volFieldCageCross">
  <materialref ref="LAr"/> 
  <solidref ref="FieldCageCross"/>
  <physvol> 
   <volumeref ref="volFieldCageBarAngle1"/>
   <position name="posFieldCageBarAngle1" unit="cm" x="0" y="0" z="0"/>
   <rotation name="rMinus45AboutZ" unit="deg" x="0" y="0" z="-45"/>
  </physvol> 
  <physvol>
   <volumeref ref="volFieldCageBarAngle2"/>
   <position name="posFieldCageBarAngle2" unit="cm" x="0" y="0" z="0"/>
   <rotation name="rPlus45AboutZ" unit="deg" x="0" y="0" z="45"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCrossBox"/>
   <position name="posFCECrossBox" unit="cm" x="0" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTriangle1"/>
   <position name="posTriangle1" unit="cm" x="(kFieldCageBeamYInt+0.5*kFieldCageCrossWidth)" y="(kFieldCageBeamYInt+0.5*kFieldCageBeamWidth)" z="0"/>
   <rotation name="rPlus0AboutZ" unit="deg" x="0" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTriangle2"/>
   <position name="posTriangle2" unit="cm" x="(kFieldCageBeamYInt+0.5*kFieldCageCrossWidth)" y="-(kFieldCageBeamYInt+0.5*kFieldCageBeamWidth)" z="0"/>
   <rotation name="rMinus90AboutZ" unit="deg" x="0" y="0" z="180"/>
  </physvol>
  <physvol>
   <volumeref ref="volTriangle1"/>
   <position name="posTriangle3" unit="cm" x="-(kFieldCageBeamYInt+0.5*kFieldCageCrossWidth)" y="-(kFieldCageBeamYInt+0.5*kFieldCageBeamWidth)" z="0"/>
   <rotation name="rPlus180AboutZ" unit="deg" x="0" y="0" z="180"/>
  </physvol>
  <physvol>
   <volumeref ref="volTriangle2"/>
   <position name="posTriangle4" unit="cm" x="-(kFieldCageBeamYInt+0.5*kFieldCageCrossWidth)" y="(kFieldCageBeamYInt+0.5*kFieldCageBeamWidth)" z="0"/>
  </physvol>
 </volume>
EOF

# Print field cage beam-loop spacers
print FIELDCAGE <<EOF;
 <volume name="volFieldCageBarFlat">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="FieldCageBarFlat"/>  
 </volume>  
 <volume name="volFieldCageBarFlat2">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="FieldCageBarFlat2"/>  
 </volume>  
 <volume name="volFieldCageBarLoopSpacer">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="FieldCageBarLoopSpacer"/>  
 </volume>  
 <volume name="volFieldCageBarLoopSpacers">
  <materialref ref="LAr"/> 
  <solidref ref="FieldCageBarLoopSpacers"/> 
EOF
if ( $spacers_on_off eq "on" ) { 
	$space=0;
	$i=1;
	while ( $space < ( $field_cage_width / 2 ) ) {
		print FIELDCAGE <<EOF;
  <physvol>
   <volumeref ref="volFieldCageBarLoopSpacer"/>
   <position name="posFieldCageBarLoopSpacerPlus$i" unit="cm" x="$space" y="0" z="0"/>
  </physvol>
EOF
		if ( $space != 0 ) {
			print FIELDCAGE <<EOF;
  <physvol>
   <volumeref ref="volFieldCageBarLoopSpacer"/>
   <position name="posFieldCageBarLoopSpacerMinus$i" unit="cm" x="-$space" y="0" z="0"/>
  </physvol>
EOF
		}
		$space+=4;
		$i++;
	}
}
print FIELDCAGE <<EOF;
 </volume>
EOF

	# Prints all components into volFieldCage
	print FIELDCAGE <<EOF;
 <volume name="volFieldCage">
  <materialref ref="LAr"/> 
  <solidref ref="FieldCage"/> 
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross1" unit="cm" x="0" y="0" z="kFieldCageCrossZPos"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross3" unit="cm" x="0" y="kFieldCageCrossYPos" z="0.8*kFieldCageBeamZInt"/>
   <rotation name="rFieldCageCross3" unit="cm" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross7" unit="cm" x="0" y="kFieldCageCrossYPos" z="-0.8*kFieldCageBeamZInt"/>
   <rotation name="rFieldCageCross3" unit="cm" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross8" unit="cm" x="0" y="0" z="kFieldCageCrossZPos"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross9" unit="cm" x="0" y="-kFieldCageCrossYPos" z="0.8*kFieldCageBeamZInt"/>
   <rotation name="rFieldCageCross9" unit="cm" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross10" unit="cm" x="0" y="-kFieldCageCrossYPos" z="0.4*kFieldCageBeamZInt"/>
   <rotation name="rFieldCageCross10" unit="cm" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross11" unit="cm" x="0" y="-kFieldCageCrossYPos" z="0*kFieldCageBeamZInt"/>
   <rotation name="rFieldCageCross11" unit="cm" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
  <position name="posFieldCageCross12" unit="cm" x="0" y="-kFieldCageCrossYPos" z="-0.4*kFieldCageBeamZInt"/>
   <rotation name="rFieldCageCross12" unit="cm" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross13" unit="cm" x="0" y="-kFieldCageCrossYPos" z="-0.8*kFieldCageBeamZInt"/>
   <rotation name="rFieldCageCross13" unit="cm" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageBarFlat"/>
   <position name="posFieldCageBar1" unit="cm" x="0" y="-kFieldCageBeamYInt" z="kFieldCageBeamZPos"/>
   <rotation name="rFieldCageBar1" unit="cm" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageBarFlat"/>
   <position name="posFieldCageBar2" unit="cm" x="0" y="0" z="kFieldCageBeamZPos"/>
   <rotation name="rFieldCageBar2" unit="cm" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageBarFlat"/>
   <position name="posFieldCageBar3" unit="cm" x="0" y="kFieldCageBeamYInt" z="kFieldCageBeamZPos"/>
   <rotation name="rFieldCageBar3" unit="cm" x="-90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageCross"/>
   <position name="posFieldCageCross2" unit="cm" x="0" y="0" z="-kFieldCageCrossZPos"/>
   <rotation name="rFieldCageCross2" unit="cm" x="0" y="180" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageBarFlat"/>
   <position name="posFieldCageBar4" unit="cm" x="0" y="-kFieldCageBeamYInt" z="-kFieldCageBeamZPos"/>
   <rotation name="rFieldCageBar4" unit="cm" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageBarFlat"/>
   <position name="posFieldCageBar5" unit="cm" x="0" y="0" z="-kFieldCageBeamZPos"/>
   <rotation name="rFieldCageBar5" unit="cm" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageBarFlat"/>
   <position name="posFieldCageBar6" unit="cm" x="0" y="kFieldCageBeamYInt" z="-kFieldCageBeamZPos"/>
   <rotation name="rFieldCageBar6" unit="cm" x="90" y="0" z="0"/>
  </physvol>
EOF

$count=7;
for ( $i=-5; $i < 6; $i++ ) {
    $fraction=$i/5;
    print FIELDCAGE <<EOF;
   <physvol>
    <volumeref ref="volFieldCageBarFlat"/>
    <position name="posFieldCageBarTop$count" unit="cm" x="0" y="kFieldCageBeamYPos+0.5*kFieldCageBeamDepth" z="$fraction*kFieldCageBeamZInt"/>
    <rotation name="rFieldCageBarTop$count" unit="cm" x="0" y="0" z="0"/>
   </physvol>
   <physvol>
    <volumeref ref="volFieldCageBarFlat2"/>
    <position name="posFieldCageBar2Top$count" unit="cm" x="0" y="(kFieldCageBeamYPos-0.25*kFieldCageBeamDepth-2*kFieldCageTubeRadius)" z="$fraction*kFieldCageBeamZInt"/>
    <rotation name="rFieldCageBar2Top$count" unit="cm" x="0" y="0" z="0"/>
   </physvol>
   <physvol>
    <volumeref ref="volFieldCageBarFlat"/>
    <position name="posFieldCageBarBot$count" unit="cm" x="0" y="-kFieldCageBeamYPos-0.5*kFieldCageBeamDepth" z="$fraction*kFieldCageBeamZInt"/>
    <rotation name="rFieldCageBarBot$count" unit="cm" x="-180" y="0" z="0"/>
   </physvol>
   <physvol>
    <volumeref ref="volFieldCageBarFlat2"/>
    <position name="posFieldCageBar2Bot$count" unit="cm" x="0" y="-(kFieldCageBeamYPos-0.25*kFieldCageBeamDepth-2*kFieldCageTubeRadius)" z="$fraction*kFieldCageBeamZInt"/>
    <rotation name="rFieldCageBar2Bot$count" unit="cm" x="-180" y="0" z="0"/>
   </physvol>
   <physvol>
    <volumeref ref="volFieldCageBarLoopSpacers"/>
    <position name="posFieldCageSpacersTop$count" unit="cm" x="0" y="0.5*kFieldCageLoopHeight-kFieldCageTubeRadius" z="$fraction*kFieldCageBeamZInt"/>
   </physvol>
   <physvol>
    <volumeref ref="volFieldCageBarLoopSpacers"/>
    <position name="posFieldCageSpacersBot$count" unit="cm" x="0" y="-(0.5*kFieldCageLoopHeight-kFieldCageTubeRadius)" z="$fraction*kFieldCageBeamZInt"/>
   </physvol>
EOF
     $count++;
}

$space=0;
$i=1;
while ( $space < ( $field_cage_width / 2 ) ) {
	$xPos=$space+2;
	print FIELDCAGE <<EOF;
  <physvol>
   <volumeref ref="volFieldCageLoop"/>
   <position name="posFieldCageLoopPlus$i" unit="cm" x="$xPos" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageLoop"/>
   <position name="posFCLoopMinus$i" unit="cm" x="-$xPos" y="0" z="0"/>
  </physvol>
EOF
	$space+=4*$field_cage_loop_interval;
	$i++;
}

print FIELDCAGE <<EOF;
 </volume>
</structure>
EOF

close(FIELDCAGE); 
