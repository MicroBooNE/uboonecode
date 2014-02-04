#!/usr/bin/perl

# This program creates GDML sub-files, with values supplied by user
# parameters.  Geometry/gdml/make_gdml.pl "zips" together those
# sub-files to make a single detector description.

# Packages
use Math::Trig;
use Math::BigFloat;
Math::BigFloat->precision(-10);
use XML::LibXML;
use Getopt::Long;

# Get the input parameters from an XML file. Optionally append a
# suffix to the GDML sub-files we create.

GetOptions( "input|i:s" => \$input,
	    "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output,
	    "wires|w:s" => \$wires);

if ( defined $help )
{
    # If the user requested help, print the usage notes and exit.
    usage();
    exit;
}

if ( ! defined $suffix )
{
    # The user didn't supply a suffix, so append nothing to the file
    # names.
    $suffix = "";
}
else
{
    # Otherwise, stick a "-" before the suffix, so that a suffix of
    # "test" applied to filename.gdml becomes "filename-test.gdml".
    $suffix = "-" . $suffix;
}


# Create an XML parser.
$parser = new XML::LibXML;

# Read in the parameters from an XML file. The following command
# slurps the entire file into a DOM data structure.
$xmldata = $parser->parse_file($input);

# Go through each parameter in the DOM data structure:
foreach $parameter ( $xmldata->findnodes('/parameters/geometry/parameter') )
{
    # Get the name and value attributes for that parameter:
    $name = $parameter->getAttribute("name");
    $value = $parameter->getAttribute("value");

    # Here's the clever part: The following eval creates a variable
    # with the same name as $name. For example, if $name eq "TPCDepth",
    # then the following statement assigns the value to $TPCDepth. The
    # value is in quotes, because some of the parameters have text
    # strings in them (like "$inch").

    eval "\$$name = '$value'";
}

# Our calculations and constants depend on the geometry of the wires.
my $SinUVAngle = sin( deg2rad($UVAngle) );
my $CosUVAngle = cos( deg2rad($UVAngle) );
my $TanUVAngle = tan( deg2rad($UVAngle) );

my $inch=2.54;

my $wires_on=1; 			# turn wires on=1 or off=0
if ( defined $wires )
{
#The user supplied the wires on parameter, so using that. Otherwise the default=1 is used.
$wires_on = $wires;
}




my $WireInterval=10;
my $NumberOfTPCPlanes=3;
my $pmt_switch="on";		#turn on or off depending on pmts wanted
my $NumberOfTestBoxes=30;
my $granite_block="off";
my $enclosureExtras="on";       #turn on or off depending on whether you'd like to generate the external things around the cryostat (ie. insulation, platform, stands, etc.) in the gdml file
my $vetoWall_switch="off";  #turn on or off a proposed scintillator wall infront of the cryostat


# The routines that create the GDML sub-files. Most of the explanatory
# comments are in gen_defs().
gen_defs();
gen_rotations();
gen_materials();

gen_microplane();
gen_microvertplane();

 gen_groundplate();	# physical volumes defined in gen_tpc()
 gen_cathode();		# physical volumes defined in gen_tpc()
 gen_fieldcage();	# physical volumes defined in gen_tpc()
gen_tpc();

if ( $pmt_switch eq "on" ) {  gen_pmt();	}	# physical volumes defined in gen_cryostat()
if ( $granite_block eq "on" ) {  gen_granite(); } # physical volumes defined in gen_enclosure()
#gen_testbox();
if ( $enclosureExtras eq "on" ) {  gen_enclosureExtras(); } #generation of insulation, etc. will happen if specified
gen_cryostat();
if ( $vetoWall_switch eq "on" ) {  gen_vetoWall();  } # physical volumes defined in gen_vetoWall()

gen_enclosure();
gen_world();
write_fragments();

exit;



sub usage()
{
    print "Usage: $0 [-h|--help] -i|--input <parameters-file> [-o|--output <fragments-file>] [-s|--suffix <string>]\n";
    print "       -i/--input can be omitted; <parameters-file> contains geometry and material parameters\n";
    print "       if -o is omitted, output goes to STDOUT; <fragments-file> is input to make_gdml.pl\n";
    print "       -s <string> appends the string to the file names; useful for multiple detector versions\n";
    print "       -h prints this message, then quits\n";
}



# Create the detector constant file. This file is actually temporary,
# since the make_gdml.pl program will interpret its contents rather
# than include it in the final GDML file.

sub gen_defs()
{
    #TPCWirePlaneWidth is the size in the y direction
    #TPCWirePlaneLength is the size in the z direction
    $TPCWirePlaneWidth	=	233;
    $TPCWirePlaneLength	=	1037;

    $pi   = pi;
    $inch = 2.54;

    $WorldWidth		=	100.0*$DetEnclosureWidth;
    $WorldHeight	=	100.0*$DetEnclosureHeight;
    $WorldLength	=	100.0*$DetEnclosureLength;

    $CathodePlateDepth	=	0.1;
    $CathodeWidth	=	5.1;
    $CathodeHeight	=	240;
    $CathodeLength	=	1042;

    $GroundPlateWidth	=	224;
    $GroundPlateHeight	=	0.1;
    $GroundPlateLength	=	1100;

    $tpc_neg_length	=	1000;
    $tpc_neg_width	=	200;
    $tpc_neg_height	=	200;
    $wire_frame_width	=	9.5;
    $wire_plane_height	=	240;
    $wire_plane_length	=	1042;
    $wires_plength		=	($wire_plane_length - 2*$wire_frame_width) ;
    $wires_pwidth		=	($wire_plane_height - 2*$wire_frame_width) ;

    $field_cage_width			=	256;
    $field_cage_height			=	207.6;  #increased size to avoid overlap of TubeBottom and Top with TPCActive
    $field_cage_cross_length	=	sqrt(($field_cage_width)**2+($field_cage_height-50)**2);
    $field_cage_length			=	1011.6; #increased size to avoid overlap of TubeFront and Back with TPCActive
    $field_cage_loop_interval	=	1; 	# =1 is normal, =4 skips 3/4

    $FieldCageTPCClearance	=	5*$inch;
    $FieldCageTubeRadius	=	0.5*$inch;
    $FieldCageTubeThickness	=	0.25*$inch;
    $FieldCageBeamDepth		=	12.5;
    $FieldCageBeamWidth		=	2.5;
    $FieldCageCrossDepth	=	2.5;
    $FieldCageCrossWidth	=	4;
    $FieldCageCrossLength	=	$field_cage_cross_length;

    $TPCTotalLength	=	$field_cage_length;
    $TPCTotalWidth	=	$field_cage_width;
    $TPCTotalHeight	=	$field_cage_height;

    $FieldCageLoopLength	=	$TPCTotalLength+2*($FieldCageTPCClearance+2*$FieldCageTubeRadius);
    $FieldCageLoopWidth		=	$field_cage_width;
    $FieldCageLoopHeight	=	$TPCTotalHeight+2*($FieldCageTPCClearance+2*$FieldCageTubeRadius);

    $FieldCageCornerRadius	=	0.5*$FieldCageTPCClearance;
    $FieldCageCornerY		=	(0.5*$FieldCageLoopHeight)-$FieldCageCornerRadius-$FieldCageTubeRadius;
    $FieldCageCornerZ		=	(0.5*$FieldCageLoopLength)-$FieldCageCornerRadius-$FieldCageTubeRadius;

    $FieldCageHeight	=	$FieldCageLoopHeight+2*($FieldCageBeamDepth+$FieldCageCrossDepth);
    $FieldCageLength	=	$FieldCageLoopLength+2*($FieldCageBeamDepth+$FieldCageCrossDepth);
    $FieldCageWidth	=	$field_cage_width;

}


sub gen_rotations()
{
    my $WirePlusRotation = $UVAngle + 90;
    my $WireMinusRotation = $UVAngle - 90;

    $ROTATIONS = "micro-rotations" . $suffix . ".gdml";
    push (@gdmlFiles, $ROTATIONS); # Add file to list of GDML fragments
    $ROTATIONS = ">" . $ROTATIONS;
    open(ROTATIONS) or die("Could not open file $ROTATIONS for writing");

    print ROTATIONS <<EOF;
<?xml version='1.0'?>
<define>
   <rotation name="rPlus30AboutX"  unit="deg" x="30"  y="0"   z="0"/>
   <rotation name="rPlus60AboutX"  unit="deg" x="60"  y="0"   z="0"/>
   <rotation name="rPlus90AboutX"  unit="deg" x="90"  y="0"   z="0"/>
   <rotation name="rMinus90AboutX"  unit="deg" x="-90"  y="0"   z="0"/>
   <rotation name="rPlusUVAngleAboutX"  unit="deg" x="90+$UVAngle" y="0"   z="0"/>
   <rotation name="rPlus150AboutX"      unit="deg" x="150" y="0"   z="0"/>
   <rotation name="rPlus180AboutX"      unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rMinusUVAngleAboutX" unit="deg" x="-30" y="0"   z="0"/>
   <rotation name="rPlus30AboutY"  unit="deg" x="0"   y="30"  z="0"/>
   <rotation name="rPlus60AboutY"  unit="deg" x="0"   y="60"  z="0"/>
   <rotation name="rPlus90AboutY"  unit="deg" x="0"   y="90"  z="0"/>
   <rotation name="rPlus180AboutY" unit="deg" x="0"   y="180" z="0"/>
   <rotation name="rMinus90AboutY" unit="deg" x="0"   y="-90" z="0"/>
   <rotation name="rPlus90AboutZ"  unit="deg" x="0"   y="0"   z="90"/>
   <rotation name="rMinus90AboutZ"  unit="deg" x="0"   y="0"   z="-90"/>
   <rotation name="rPlus180AboutZ"      unit="deg" x="0"   y="0"   z="180"/>
   <rotation name="rMinus180AboutZ"     unit="deg" x="0"   y="0"   z="-180"/>
   <rotation name="rMinus90AboutYPlus180AboutZ" unit="deg" x="0" y="-90" z="180"/>
   <rotation name="rMinus90AboutYMinus90AboutZ" unit="deg" x="0" y="-90" z="-90"/>
   <rotation name="rPlus90AboutYPlus180AboutZ" unit="deg" x="0" y="90" z="180"/>
   <rotation name="rMinus90AboutYPlus90AboutZ" unit="deg" x="0" y="-90" z="90"/>
   <rotation name="rPlus90AboutYMinus90AboutZ" unit="deg" x="0" y="90" z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutZ"  unit="deg" x="90" y="0"   z="90"/>
   <rotation name="rPlus90AboutXPlus180AboutZ" unit="deg" x="90" y="0"   z="180"/>
   <rotation name="rPlus90AboutXMinus90AboutY" unit="deg" x="90" y="-90" z="0"/>
   <rotation name="rPlus90AboutXMinus90AboutZ" unit="deg" x="90" y="0"   z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutY"  unit="deg"  x="90" y="90" z="0"/>
   <rotation name="rPMTRotation1"  unit="deg" x="90"  y="270"   z="0"/>
   <position name="posCenter" unit="mm" x="0" y="0" z="0"/>

</define>
EOF
    close (ROTATIONS);
}


sub gen_materials()
{
    # Create the materials file name and open it.
    $MATERIALS = "materials" . $suffix . ".gdml";
    push (@gdmlFiles, $MATERIALS); # Add file to list of GDML fragments
    $MATERIALS = ">" . $MATERIALS;
    open(MATERIALS) or die("Could not open file $MATERIALS for writing");

    # Write the standard XML prefix.
    print MATERIALS <<EOF;
<?xml version='1.0'?>
EOF

    # Go back the DOM structure read in near the beginning of the
    # program. For each <materials /> element (and there'll probably
    # be only one):
    foreach $materials ( $xmldata->findnodes('/parameters/materials') )
    {
	# Convert that element back to text, and write it out.
	print MATERIALS $materials->toString;
    }

    close (MATERIALS);
}


# This is a re-write of Brian Rebel's gen_microvertplane.C into
# Perl. It contructs the TPC wire plane for the Y view.

sub gen_microvertplane()
{
  my $NumberWires = 0;
  
  if ( $wires_on == 1 ) 
    { $NumberWires = int( $TPCWirePlaneLength / $TPCWirePitch ); }
  
    

    $GDML = "micro-vertplane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Define the solids and structures: the wires and the TPC wire plane.
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
<tube name="TPCWireVert"
  rmax="0.5*$TPCWireThickness"
  z="$TPCWirePlaneWidth"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="TPCPlaneVert"
  x="$TPCWirePlaneThickness" 
  y="$TPCWirePlaneWidth" 
  z="$TPCWirePlaneLength"
  lunit="cm"/>
</solids>
<structure>
  <volume name="volTPCWireVert">
    <materialref ref="Titanium"/>
    <solidref ref="TPCWireVert"/>
  </volume>
  <volume name="volTPCPlaneVert">
    <materialref ref="LAr"/>       
    <solidref ref="TPCPlaneVert"/>
EOF

    # the wires 
    for ( $i = 0; $i < $NumberWires; ++$i)
    {
	print GDML <<EOF;
	    <physvol>
	     <volumeref ref="volTPCWireVert"/>
	     <position name="posTPCWireVert$i" unit="cm" z="-0.5*$TPCWirePlaneLength+$TPCWirePitch*($i+1)" x="0" y="0"/>
	     <rotationref ref="rPlus90AboutX"/>
	    </physvol>
EOF
    }

    print GDML <<EOF;
  </volume>
</structure>
</gdml>
EOF

    close(GDML);
}


# This is a re-write of Brian Rebel's gen_microplane.C into Perl. It
# constructs the TPC wire plane for the U or V view.

sub gen_microplane()
{

    $TPCWireXPitch=$TPCWirePitch/$CosUVAngle;

    $GDML = "micro-plane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Calculate the number of wire ends on a given y-edge of the plane.
    my $TPCYWirePitch = $TPCWirePitch / $CosUVAngle;

   my $NumberWiresPerEdge = 0;
   if ( $wires_on == 1 ) 
    {  $NumberWiresPerEdge = int( $TPCWirePlaneLength / $TPCYWirePitch );}
  
    # How many side wires will be "cut off" by the lower or higher
    # z-edge?
   my $NumberSideWires = 0;
    if ( $wires_on == 1 ) 
    { $NumberSideWires = int( $TanUVAngle * $TPCWirePlaneWidth / $TPCYWirePitch ); }

    # The number of full-length "center" wires.
  
   my $NumberCenterWires = 0.;
   if ( $wires_on == 1 ) 
    { $NumberCenterWires = $NumberWiresPerEdge - $NumberSideWires; }

    # define the solids
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
EOF

    # wires on either end of the tpc #subtraction added to avoid overlap for all wires
    for($i = 0; $i < $NumberSideWires; ++$i)
    {
	print GDML <<EOF;
<tube name="TPCWire$i"
  rmax="0.5*$TPCWireThickness"
  z="$TPCWireXPitch*($i+1)/$SinUVAngle-0.01498648"    
  deltaphi="360"
  aunit="deg"
  lunit="cm"/> 
EOF
    }

    # The solids for the middle wire and the TPC wire plane, and start off the structures.
    #subtraction added to avoid overlap for all wires
    print GDML <<EOF;
<tube name="TPCWireCommon"
  rmax="0.5*$TPCWireThickness"
  z="$TPCWirePlaneWidth/$CosUVAngle-0.0259573526"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="TPCPlane"
  x="$TPCWirePlaneThickness"
  y="$TPCWirePlaneWidth"
  z="$TPCWirePlaneLength"
  lunit="cm"/>
</solids>
<structure>
EOF
 
    # the wires at either end of the plane
    for ($i = 0; $i < $NumberSideWires; ++$i)
    {
	print GDML <<EOF;
    <volume name="volTPCWire$i">
	<materialref ref="Titanium"/>
	<solidref ref="TPCWire$i"/>
    </volume>
EOF
    }

  
    # The wires in the middle of the plane, and the plane itself.
    print GDML <<EOF;
    <volume name="volTPCWireCommon">
      <materialref ref="Titanium"/>
      <solidref ref="TPCWireCommon"/>
    </volume>
    <volume name="volTPCPlane">
      <materialref ref="LAr"/>
      <solidref ref="TPCPlane"/>
EOF

  # the wires at the -z end
  for ($i = 0; $i < $NumberSideWires; ++$i)
  {
    $j=$NumberSideWires-$i-1;

    print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWire$j"/> 
     <position name="posTPCWireF$j" unit="cm" y="-0.5*($i)*$TPCWireXPitch/$TanUVAngle" z="-0.5*$TPCWirePlaneLength+0.5*$TPCWireXPitch*($j+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  $ypos=-0.5*($i)*$TPCWireXPitch/$TanUVAngle;
  $zpos=-0.5*$TPCWirePlaneLength+0.5*$TPCWireXPitch*($j);
  open (MYFILE, '>>data.txt');
  print MYFILE "TPCWire$j y=$ypos z=$zpos\n";
  }

  # the wires at the +z end
  for ($i = 0; $i < $NumberSideWires; ++$i)
  {
      my $j = $NumberSideWires-$i-1;
      my $k = $NumberCenterWires+$NumberSideWires+$i;
      $zposlast=-0.5*$TPCWirePlaneLength+$TPCWireXPitch*(0.5*$NumberSideWires+$NumberCenterWires);
      $zpos=$zposlast+0.5*$TPCWireXPitch*($i);

    print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWire$j"/> 
     <position name="posTPCWireB$j" unit="cm" y="0.5*($i)*$TPCWireXPitch/$TanUVAngle" z="$zpos" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF


  }

  # The wires in the middle.
  for ($i = 0; $i < $NumberCenterWires - 1; ++$i)
  {
      my $j = $NumberSideWires+$i;
      $zpos=-0.5*$TPCWirePlaneLength+$TPCWireXPitch*(0.5*$NumberSideWires+$i+1);

      print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWireCommon"/>
     <position name="posTPCWire$j" unit="cm" y="0" z="$zpos" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  }


      print GDML <<EOF;
  </volume>
</structure>
</gdml>
EOF

  close(GDML);
}


# subdirectory to write field cage
sub gen_fieldcage() {

    # Set up the output file.
    $FIELDCAGE = "micro-fieldcage.gdml";
    push (@gdmlFiles, $FIELDCAGE); # Add file to list of constant files
    $FIELDCAGE = ">" . $FIELDCAGE;
    open(FIELDCAGE) or die("Could not open file $FIELDCAGE for writing");

    # set the field cage constants

    # Prints Field Cage solids
    print FIELDCAGE <<EOF;
<solids>
 <tube name="FieldCageTubeZ"
  rmin="$FieldCageTubeRadius-$FieldCageTubeThickness"
  rmax="$FieldCageTubeRadius"  
  z="$FieldCageLoopLength-2*($FieldCageCornerRadius+$FieldCageTubeRadius)" 
  deltaphi="360" 
  aunit="deg" 
  lunit="cm"/> 
 <tube name="FieldCageTubeY" 
  rmin="$FieldCageTubeRadius-$FieldCageTubeThickness"
  rmax="$FieldCageTubeRadius"  
  z="$FieldCageLoopHeight-2*($FieldCageCornerRadius+$FieldCageTubeRadius)"  
  deltaphi="360" 
  aunit="deg" 
  lunit="cm"/> 
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
</structure>
EOF

close(FIELDCAGE); 

}


# Cathode solids and volumes
sub gen_cathode() {

    $CATHODE = "micro-cathode" . $suffix . ".gdml";
    push (@gdmlFiles, $CATHODE); # Add file to list of GDML fragments
    $CATHODE = ">" . $CATHODE;
    open(CATHODE) or die("Could not open file $CATHODE for writing");

    print CATHODE <<EOF;

<?xml version='1.0'?>
<gdml>
<solids>
 <box name="CathodePlate"
  lunit="cm"
  x="$CathodePlateDepth"
  y="$CathodeHeight"
  z="$CathodeLength"/>
</solids>
<structure>
 <volume name="volCathodePlate">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="CathodePlate"/>
 </volume>
</structure>
</gdml>
EOF

}

sub gen_groundplate() {
#this subroutine will produce the gdml fragment for ground plate

    $GROUNDPLATE = "micro-groundplate" . $suffix . ".gdml";
    push (@gdmlFiles, $GROUNDPLATE); # Add file to list of GDML fragments
    $GROUNDPLATE = ">" . $GROUNDPLATE;
    open(GROUNDPLATE) or die("Could not open file $GROUNDPLATE for writing");

    print GROUNDPLATE <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="GroundPlate"
  lunit="cm"
  x="$GroundPlateWidth"
  y="$GroundPlateHeight"
  z="$GroundPlateLength"/>
</solids>
<structure>
 <volume name="volGroundPlate">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="GroundPlate"/>
 </volume>
</structure>
</gdml>
EOF

    close(GROUNDPLATE);
}


# Parameterize the TPC and the planes within it.
sub gen_tpc()
{
    # Set up the output file.
    $GDML = "micro-tpc" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");


    # Size info for active TPC volume (LAr inside)
    $aTPC_xos_cathode = 0.5*$CathodeWidth + 0.5*$CathodePlateDepth;
    $aTPC_xos_wires = 3*$TPCWirePlaneSpacing + $TPCWirePlaneThickness;
    $aTPC_xoffset = $aTPC_xos_cathode - $aTPC_xos_wires ; 
    $TPCActiveDepth  = $TPCDepth - $aTPC_xos_cathode - $aTPC_xos_wires ; 
#    $TPCActiveHeight = $FieldCageLoopHeight - 4*$FieldCageTubeRadius ;   
#    $TPCActiveLength = $FieldCageLoopLength - 4*$FieldCageTubeRadius ;   
    $TPCActiveHeight = $TPCWirePlaneWidth;   #  
    $TPCActiveLength = $TPCWirePlaneLength-0.2;   # extra subtraction to arrive at TPCActive values in the TDR


    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="TPC"
  lunit="cm"
  x="$TPCDepth"
  y="$TPCWidth"
  z="$TPCLength"/>
 <box name="TPCActive"
  lunit="cm"
  x="$TPCActiveDepth"
  y="$TPCActiveHeight"
  z="$TPCActiveLength"/>
</solids>
<structure>
 <volume name="volTPCActive">
  <materialref ref="LAr"/>
  <solidref ref="TPCActive"/>
 </volume>
 <volume name="volTPC">
  <materialref ref="LAr"/>
  <solidref ref="TPC"/>
EOF

    # Ground Plane physical volumes
    # Center = (9,124,0)
    $ground_plate_X=9;
    $ground_plate_Y=124;
    $ground_plate_Z=0;

    print GDML <<EOF;
<!--  <physvol>
   <volumeref ref="volGroundPlate"/>
   <position name="posGroundPlate" unit="cm" x="$ground_plate_X+0.25" y="$ground_plate_Y" z="0"/>
  </physvol>-->
EOF


    # Cathode Plane physical volumes
    # Center = (0.5*($TPCDepth-$CathodeWidth),0,0)
    $Cathode_X=0.5*($TPCDepth-$CathodeWidth);
    $Cathode_plate_offset=0.5*$CathodePlateDepth;
    $Cathode_frame_position=$Cathode_X+$Cathode_plate_offset;

    print GDML <<EOF;
  <physvol>
   <volumeref ref="volCathodePlate"/>
   <position name="posCathodePlate" unit="cm" x="$Cathode_X" y="0" z="0"/>
  </physvol>
EOF


    # Wire Plane physical volumes 
    # Center = ( -0.5*($TPCWidth-kWireFrameWidth) , 0 , 0 )
    print GDML <<EOF;
  <physvol>
   <volumeref ref="volTPCPlaneVert"/>
   <position name="posTPCPlaneVert" unit="cm" x="-0.5*$TPCActiveDepth-2*$TPCWirePlaneSpacing" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCPlane"/>
   <position name="posTPCPlane" unit="cm" x="-0.5*$TPCActiveDepth-$TPCWirePlaneSpacing" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCPlane"/>
   <position name="posTPCPlane2" unit="cm" x="-0.5*$TPCActiveDepth" y="0" z="0"/>
   <rotationref ref="rPlus180AboutY"/>
  </physvol>
EOF


$space=0;
$i=1;
while ( $space < ( $field_cage_width / 2 ) ) {
	$xPos=$space+2;
	print GDML <<EOF;
  <physvol>
   <volumeref ref="volFieldCageTubeTop"/>
   <position name="posFieldCageTubeTopA$i" unit="cm" x="$xPos" y="0.5*($FieldCageLoopHeight-2*$FieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBot"/>
   <position name="posFieldCageTubeBotA$i" unit="cm" x="$xPos" y="-0.5*($FieldCageLoopHeight-2*$FieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeFront"/>
   <position name="posFieldCageTubeFrontA$i" unit="cm" x="$xPos" y="0" z="0.5*($FieldCageLoopLength-2*$FieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertPlusA$i" unit="deg" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBack"/>
   <position name="posFieldCageTubeBackA$i" unit="cm" x="$xPos" y="0" z="-0.5*($FieldCageLoopLength-2*$FieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertMinusA$i" unit="deg" x="-90" y="0" z="0"/>
  </physvol>
EOF
	print GDML <<EOF;
  <physvol>
   <volumeref ref="volFieldCageTubeTop"/>
   <position name="posFieldCageTubeTopB$i" unit="cm" x="-$xPos" y="0.5*($FieldCageLoopHeight-2*$FieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBot"/>
   <position name="posFieldCageTubeBotB$i" unit="cm" x="-$xPos" y="-0.5*($FieldCageLoopHeight-2*$FieldCageTubeRadius)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeFront"/>
   <position name="posFieldCageTubeFrontB$i" unit="cm" x="-$xPos" y="0" z="0.5*($FieldCageLoopLength-2*$FieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertFrontB$i" unit="deg" x="90" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volFieldCageTubeBack"/>
   <position name="posFieldCageTubeBackB$i" unit="cm" x="-$xPos" y="0" z="-0.5*($FieldCageLoopLength-2*$FieldCageTubeRadius)"/>
   <rotation name="rFieldCageVertBackB$i" unit="deg" x="-90" y="0" z="0"/>
  </physvol>
EOF
	$space+=4*$field_cage_loop_interval;
	$i++;
}


    print GDML <<EOF;
   <physvol>
    <volumeref ref="volTPCActive"/>
    <position name="posTPCActive" unit="cm" x="-$aTPC_xoffset" y="0" z="0"/>
   </physvol>
EOF

    # Closes TPC volume definition space
    print GDML <<EOF;
 </volume> 
</structure>
</gdml>
EOF

   close(GDML);
}


# Generates Ben Jones's PMT micro-pmtdef (with temporary edit to ellipsoid shapes
sub gen_pmt {

    $PMT = "micro-pmtdef" . $suffix . ".gdml";
    push (@gdmlFiles, $PMT); # Add file to list of GDML fragments
    $PMT = ">" . $PMT;
    open(PMT) or die("Could not open file $PMT for writing");

	print PMT <<EOF;
<solids>
 <tube name="PMTVolume"
  rmax="(6.1*2.54)"
  z="(11.1*2.54)+0.005"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>

 <tube name="PMT_AcrylicPlate"
  rmax="(6.0*2.54)"
  z="(0.2)"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_Stalk"
  rmax="(1.25*2.54)"
  z="(3.0*2.54)"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_SteelBase"
  rmax="(6.0*2.54)"
  z="(1.5*2.54)"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
 <tube name="PMT_Underside"
  rmax="2.54*4.0"
  z="2.54*2.5"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
EOF
	print PMT <<EOF;
 <tube name="PMT_Lens"
  rmax="2.54*4.0"
  z="2.54*2.5"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
EOF

	print PMT <<EOF;
</solids>
<structure>
 <volume name="volOpDetSensitive">
  <materialref ref="LAr"/>
  <solidref ref="PMT_AcrylicPlate"/>
 </volume>
 <volume name="vol_PMT_AcrylicPlate">
  <materialref ref="Acrylic"/>
  <solidref ref="PMT_AcrylicPlate"/>
 </volume>
 <volume name="vol_PMT_Stalk">
  <materialref ref="Glass"/>
  <solidref ref="PMT_Stalk"/>
 </volume>
 <volume name="vol_PMT_SteelBase">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="PMT_SteelBase"/>
 </volume>
 <volume name="vol_PMT_Underside">
  <materialref ref="Glass"/>
  <solidref ref="PMT_Underside"/>
 </volume>
EOF
	print PMT <<EOF;
 <volume name="vol_PMT_Lens">
  <materialref ref="LAr"/>
  <solidref ref="PMT_Lens"/>
 </volume>
EOF
	print PMT <<EOF;
 <volume name="volPMT">
  <materialref ref="LAr"/>
  <solidref ref="PMTVolume"/>
  <physvol>
   <volumeref ref="volOpDetSensitive"/>
   <position name="posOpDetSensitive" unit="cm" x="0" y="0" z="(5.5 * 2.54) - 0.1"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_AcrylicPlate"/>
   <position name="pos_PMT_AcrylicPlate" unit="cm" x="0" y="0" z="(5.5 * 2.54) - 0.3"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Stalk"/>
   <position name="pos_PMT_Stalk" unit="cm" x="0" y="0" z="(3.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_SteelBase"/>
   <position name="pos_PMT_SteelBase" unit="cm" x="0" y="0" z="(0.75 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Lens"/>
   <position name="pos_PMT_Lens" unit="cm" x="0" y="0" z="(7.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
  <physvol>
   <volumeref ref="vol_PMT_Underside"/>
   <position name="pos_PMT_Underside" unit="cm" x="0" y="0" z="(7.0 * 2.54)-(5.5 * 2.54)"/>
  </physvol>
EOF

	print PMT <<EOF;
 </volume>
</structure>
EOF

}


# Generates Bill Seligman's granite block
sub gen_granite() {

    $GRANITE = "micro-graniteblock" . $suffix . ".gdml";
    push (@gdmlFiles, $GRANITE); # Add file to list of GDML fragments
    $GRANITE = ">" . $GRANITE;
    open(GRANITE) or die("Could not open file $GRANITE for writing");

	print GRANITE <<EOF;
<solids>
  <!-- 06-Dec-2011 WGS:
       The dimensions of the granite block are
       input parameters and subject to change -->
  <box name="GraniteBlock" lunit="cm"
    x="160" y="250" z="180" />
    
  <!-- 06-Dec-2011 WGS: Cut-and-paste from the GDML description of the Double Chooz Outer Veto module -->  
    <box lunit="mm" name="ScintPlastic_32250x1bf1130" x="49.6"      y="9.6"     z="3225"/>
    <box lunit="mm" name="ScintStrip_32250x1bf10e0"   x="50"        y="10"      z="3225.2"/>
    <box lunit="mm" name="TapeSheet_32250x1beee50"    x="1628.25"   y="0.9999"  z="3225.2"/>
    <box lunit="mm" name="SideSpacer_32250x1bef410"   x="10"        y="24.4064" z="3605.2"/>
    <box lunit="mm" name="AirBox0x1bede80"            x="1648.25"   y="24.4064" z="380"/>
    <box lunit="mm" name="Module_32250x1bf1000"       x="1649.0628" y="25.2192" z="3606.0128"/>
    <!-- needed for rotation -->
    <box lunit="cm" name="Module_3225_rot"            x="2.6"       y="165"     z="361"/>
</solids>
<structure>
    <!-- 06-Dec-2011 WGS:
     Cut-and-paste from the GDML description of the Double Chooz Outer Veto module,
     with a change: I assume we don't need multiple definitions of "Air" -->
    <volume name="ScintPlastic_3225_log0x1bedfc0">
      <materialref ref="plastic0x1bec080"/>
      <solidref ref="ScintPlastic_32250x1bf1130"/>
    </volume>
    <volume name="ScintStrip_3225_log0x1beded0">
      <materialref ref="TitaniumDioxide0x1beaf50"/>
      <solidref ref="ScintStrip_32250x1bf10e0"/>
      <physvol name="ScintPlastic_3225_phy0x1bee0b0">
        <volumeref ref="ScintPlastic_3225_log0x1bedfc0"/>
      </physvol>
    </volume>
    <volume name="TapeSheet_3225_log0x1beeee0">
      <materialref ref="foamtape0x1bec350"/>
      <solidref ref="TapeSheet_32250x1beee50"/>
    </volume>
    <volume name="SideSpacer_3225_log0x1bef4a0">
      <materialref ref="rubber0x1beccc0"/>
      <solidref ref="SideSpacer_32250x1bef410"/>
    </volume>
    <volume name="AirBox_log0x1bf0810">
      <materialref ref="Air"/>
      <solidref ref="AirBox0x1bede80"/>
    </volume>
    <volume name="Module_3225_log0x1bf1050">
      <materialref ref="Aluminum0x1beb1f0"/>
      <solidref ref="Module_32250x1bf1000"/>
      <physvol name="phys_strip_3225_00x1bee160">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_00x1bee160_pos" unit="mm" x="789.075" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_10x1bee1f0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_10x1bee1f0_pos" unit="mm" x="738.975" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_20x1bee280">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_20x1bee280_pos" unit="mm" x="688.875" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_30x1bee340">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_30x1bee340_pos" unit="mm" x="638.775" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_40x1bee3d0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_40x1bee3d0_pos" unit="mm" x="588.675" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_50x1bee4b0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_50x1bee4b0_pos" unit="mm" x="538.575" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_60x1bee540">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_60x1bee540_pos" unit="mm" x="488.475" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_70x1bee5d0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_70x1bee5d0_pos" unit="mm" x="438.375" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_80x1bee660">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_80x1bee660_pos" unit="mm" x="388.275" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_90x1bee420">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_90x1bee420_pos" unit="mm" x="338.175" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_100x1b49940">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_100x1b49940_pos" unit="mm" x="288.075" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_110x1b499d0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_110x1b499d0_pos" unit="mm" x="237.975" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_120x1b49a60">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_120x1b49a60_pos" unit="mm" x="187.875" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_130x1b49af0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_130x1b49af0_pos" unit="mm" x="137.775" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_140x1b49b80">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_140x1b49b80_pos" unit="mm" x="87.675" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_150x1b49c10">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_150x1b49c10_pos" unit="mm" x="37.575" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_160x1b49ca0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_160x1b49ca0_pos" unit="mm" x="-12.525" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_170x1b49870">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_170x1b49870_pos" unit="mm" x="-62.625" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_180x1b49e40">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_180x1b49e40_pos" unit="mm" x="-112.725" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_190x1b49ed0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_190x1b49ed0_pos" unit="mm" x="-162.825" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_200x1b49f60">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_200x1b49f60_pos" unit="mm" x="-212.925" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_210x1b49ff0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_210x1b49ff0_pos" unit="mm" x="-263.025" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_220x1b4a080">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_220x1b4a080_pos" unit="mm" x="-313.125" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_230x1b4a110">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_230x1b4a110_pos" unit="mm" x="-363.225" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_240x1b4a1a0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_240x1b4a1a0_pos" unit="mm" x="-413.325" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_250x1b4a230">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_250x1b4a230_pos" unit="mm" x="-463.425" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_260x1b4a2c0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_260x1b4a2c0_pos" unit="mm" x="-513.525" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_270x1b4a350">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_270x1b4a350_pos" unit="mm" x="-563.625" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_280x1b4a3e0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_280x1b4a3e0_pos" unit="mm" x="-613.725" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_290x1b4a470">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_290x1b4a470_pos" unit="mm" x="-663.825" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_300x1b4a500">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_300x1b4a500_pos" unit="mm" x="-713.925" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_310x1b4a590">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_310x1b4a590_pos" unit="mm" x="-764.025" y="-6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_320x1b4a620">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_320x1b4a620_pos" unit="mm" x="764.025" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_330x1b49d30">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_330x1b49d30_pos" unit="mm" x="713.925" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_340x1b4a880">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_340x1b4a880_pos" unit="mm" x="663.825" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_350x1b4a8d0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_350x1b4a8d0_pos" unit="mm" x="613.725" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_360x1b4a960">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_360x1b4a960_pos" unit="mm" x="563.625" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_370x1b4a9f0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_370x1b4a9f0_pos" unit="mm" x="513.525" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_380x1b4aa80">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_380x1b4aa80_pos" unit="mm" x="463.425" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_390x1b4ab10">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_390x1b4ab10_pos" unit="mm" x="413.325" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_400x1b4aba0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_400x1b4aba0_pos" unit="mm" x="363.225" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_410x1b4ac30">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_410x1b4ac30_pos" unit="mm" x="313.125" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_420x1b4acc0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_420x1b4acc0_pos" unit="mm" x="263.025" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_430x1b4ad50">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_430x1b4ad50_pos" unit="mm" x="212.925" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_440x1b4ade0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_440x1b4ade0_pos" unit="mm" x="162.825" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_450x1b4ae70">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_450x1b4ae70_pos" unit="mm" x="112.725" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_460x1b4af00">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_460x1b4af00_pos" unit="mm" x="62.625" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_470x1b4af90">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_470x1b4af90_pos" unit="mm" x="12.525" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_480x1b4b020">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_480x1b4b020_pos" unit="mm" x="-37.575" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_490x1b4b0b0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_490x1b4b0b0_pos" unit="mm" x="-87.675" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_500x1bee6b0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_500x1bee6b0_pos" unit="mm" x="-137.775" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_510x1bee740">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_510x1bee740_pos" unit="mm" x="-187.875" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_520x1bee7d0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_520x1bee7d0_pos" unit="mm" x="-237.975" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_530x1bee860">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_530x1bee860_pos" unit="mm" x="-288.075" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_540x1bee8f0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_540x1bee8f0_pos" unit="mm" x="-338.175" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_550x1bee980">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_550x1bee980_pos" unit="mm" x="-388.275" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_560x1beea10">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_560x1beea10_pos" unit="mm" x="-438.375" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_570x1beeaa0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_570x1beeaa0_pos" unit="mm" x="-488.475" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_580x1beeb30">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_580x1beeb30_pos" unit="mm" x="-538.575" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_590x1beebc0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_590x1beebc0_pos" unit="mm" x="-588.675" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_600x1beec50">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_600x1beec50_pos" unit="mm" x="-638.775" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_610x1beece0">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_610x1beece0_pos" unit="mm" x="-688.875" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_620x1beed70">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_620x1beed70_pos" unit="mm" x="-738.975" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_strip_3225_630x1beee00">
        <volumeref ref="ScintStrip_3225_log0x1beded0"/>
        <position name="phys_strip_3225_630x1beee00_pos" unit="mm" x="-789.075" y="6.2032" z="190"/>
      </physvol>
      <physvol name="phys_tapesheet_3225_00x1beefb0">
        <volumeref ref="TapeSheet_3225_log0x1beeee0"/>
        <position name="phys_tapesheet_3225_00x1beefb0_pos" unit="mm" x="0" y="11.7032" z="190"/>
      </physvol>
      <physvol name="phys_tapesheet_3225_10x1b4a6b0">
        <volumeref ref="TapeSheet_3225_log0x1beeee0"/>
        <position name="phys_tapesheet_3225_10x1b4a6b0_pos" unit="mm" x="0" y="0.7032" z="190"/>
      </physvol>
      <physvol name="phys_tapesheet_3225_20x1b4a740">
        <volumeref ref="TapeSheet_3225_log0x1beeee0"/>
        <position name="phys_tapesheet_3225_20x1b4a740_pos" unit="mm" x="0" y="-0.7032" z="190"/>
      </physvol>
      <physvol name="phys_tapesheet_3225_30x1b4a7d0">
        <volumeref ref="TapeSheet_3225_log0x1beeee0"/>
        <position name="phys_tapesheet_3225_30x1b4a7d0_pos" unit="mm" x="0" y="-11.7032" z="190"/>
      </physvol>
      <physvol name="phys_spacer3225_10x1bef570">
        <volumeref ref="SideSpacer_3225_log0x1bef4a0"/>
        <position name="phys_spacer3225_10x1bef570_pos" unit="mm" x="819.125" y="0" z="0"/>
      </physvol>
      <physvol name="phys_spacer3225_20x1bef600">
        <volumeref ref="SideSpacer_3225_log0x1bef4a0"/>
        <position name="phys_spacer3225_20x1bef600_pos" unit="mm" x="-819.125" y="0" z="0"/>
      </physvol>
      <physvol name="phys_airbox0x1bef650">
        <volumeref ref="AirBox_log0x1bf0810"/>
        <position name="phys_airbox0x1bef650_pos" unit="mm" x="0" y="0" z="-1612.6"/>
      </physvol>
    </volume>
  
    <!-- 06-Dec-11 WGS: 
      To rotate the scintillator modules properly, we need two separate
      rotations. Take this opportunity to define the scintillator modules with 
      four distinct names. We'll need to identify the volumes separately when 
      we write hit-processing code for Geant4. -->
    <volume name="volDoubleChoozScintModule_0">
      <materialref ref="Aluminum0x1beb1f0"/>
      <solidref ref="Module_3225_rot"/>
    	<physvol>
    		<volumeref ref="Module_3225_log0x1bf1050" />
            <positionref ref="posCenter" />
            <rotationref ref="rPlus90AboutZ" />
    	</physvol>
    </volume>
    <volume name="volDoubleChoozScintModule_1">
      <materialref ref="Aluminum0x1beb1f0"/>
      <solidref ref="Module_3225_rot"/>
    	<physvol>
    		<volumeref ref="Module_3225_log0x1bf1050" />
            <positionref ref="posCenter" />
            <rotationref ref="rPlus90AboutZ" />
    	</physvol>
    </volume>
    <volume name="volDoubleChoozScintModule_2">
      <materialref ref="Aluminum0x1beb1f0"/>
      <solidref ref="Module_3225_rot"/>
    	<physvol>
    		<volumeref ref="Module_3225_log0x1bf1050" />
            <positionref ref="posCenter" />
            <rotationref ref="rPlus90AboutZ" />
    	</physvol>
    </volume>
    <volume name="volDoubleChoozScintModule_3">
      <materialref ref="Aluminum0x1beb1f0"/>
      <solidref ref="Module_3225_rot"/>
    	<physvol>
    		<volumeref ref="Module_3225_log0x1bf1050" />
            <positionref ref="posCenter" />
            <rotationref ref="rPlus90AboutZ" />
    	</physvol>
    </volume>
   <volume name="volGraniteBlock">
     <materialref ref="Granite" />
     <solidref ref="GraniteBlock" />
   </volume>

</structure>
EOF
}


#Parameterize the steel cryostat that encloses the TPC.
sub gen_cryostat()
{
    # Set up the output file.
    $CRYOSTAT = "micro-cryostat" . $suffix . ".gdml";
    push (@gdmlFiles, $CRYOSTAT); # Add file to list of GDML fragments
    $CRYOSTAT = ">" . $CRYOSTAT;
    open(CRYOSTAT) or die("Could not open file $CRYOSTAT for writing");

    print CRYOSTAT <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <tube name="Cryostat" 
  rmax="$CryostatOuterRadius"
  z="$CryostatLength+2*$CryostatEndcapThickness"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<tube name="SteelTube"
  rmin="$CryostatInnerRadius"
  rmax="$CryostatOuterRadius-0.1"
  z="$CryostatLength"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<sphere name="EndCap" rmin="144*2.54" rmax="144.5*2.54" deltaphi="360" deltatheta="31.3822" aunit="deg" lunit="cm"/>
EOF


	print CRYOSTAT <<EOF;
</solids>

<structure>
<volume name="volEndCap">
   <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
   <solidref ref="EndCap"/>
  </volume>
 <volume name="volSteelTube">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="SteelTube"/>
 </volume>
 <volume name="volCryostat">
  <materialref ref="LAr"/>
  <solidref ref="Cryostat"/>
  <physvol>
   <volumeref ref="volSteelTube"/>
   <position name="posSteelTube" unit="cm" x="0" y="0" z="0"/>
  </physvol>
<physvol>
   <volumeref ref="volEndCap"/>
   <position name="posEndCap1" unit="cm" x="0" y="0" z="427.75*2.54/2- 2.54*sqrt(144.5^2-75.5^2)"/>
   </physvol>
   <physvol>
    <volumeref ref="volEndCap"/>
    <position name="posEndCap2" unit="cm" x="0" y="0" z="-(427.75*2.54/2 - 2.54*sqrt(144.5^2-75.5^2))"/>
    <rotationref ref="rPlus180AboutY"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPC"/>
   <position name="posTPC" unit="cm" x="0.0" y="0.97" z="0"/>
  </physvol>
EOF


  @pmt_pos = ( ' x="-147.8"  y="3.21654"  z="-472"',
               ' x="-147.76" y="-52.6635" z="-420"',
               ' x="-147.8"  y="59.0965"  z="-420"',
               ' x="-147.76" y="-52.6635" z="-380"',
               ' x="-147.8"  y="59.0965"  z="-380"',
               ' x="-147.8"  y="3.21654"  z="-328"',
               ' x="-147.8"  y="3.21654"  z="-272"',
               ' x="-147.76" y="-52.6635" z="-220"',
               ' x="-147.8"  y="59.0965"  z="-220"',
               ' x="-147.76" y="-52.6635" z="-180"',
               ' x="-147.8"  y="59.0965"  z="-180"',
               ' x="-147.8"  y="3.21654"  z="-128"',
               ' x="-147.8"  y="3.21654"  z="-72"',
               ' x="-147.76" y="-52.6635" z="-20"',
               ' x="-147.8"  y="59.0965"  z="-20"',
               ' x="-147.76" y="-52.6635" z="20"',
               ' x="-147.8"  y="59.0965"  z="20"',
               ' x="-147.8"  y="3.21654"  z="72"',
               ' x="-147.8"  y="3.21654"  z="128"',
               ' x="-147.76" y="-52.6635" z="180"',
               ' x="-147.8"  y="59.0965"  z="180"',
               ' x="-147.76" y="-52.6635" z="220"',
               ' x="-147.8"  y="59.0965"  z="220"',
               ' x="-147.8"  y="3.21654"  z="272"',
               ' x="-147.8"  y="3.21654"  z="328"',
               ' x="-147.76" y="-52.6635" z="380"',
               ' x="-147.8"  y="59.0965"  z="380"',
               ' x="-147.76" y="-52.6635" z="420"',
               ' x="-147.8"  y="59.0965"  z="420"',
               ' x="-147.8"  y="3.21654"  z="472"' );

  if ( $pmt_switch eq "on" ) {
    for ( $i=0; $i<30; ++$i ){
      print CRYOSTAT <<EOF;
  <physvol>
   <volumeref ref="volPMT"/>
   <position name="posPMT$i" unit="cm" @pmt_pos[$i]/>
   <rotationref ref="rPMTRotation1"/>
  </physvol>
EOF
    }
  }

  if ( $test_switch eq "on" ) {
    for ( $i=0; $i<$NumberOfTestBoxes; ++$i ){
      print CRYOSTAT <<EOF;
  <physvol>
   <volumeref ref="volTEST"/>
   <position name="posTEST$i" unit="cm" @pmt_pos[$i]/>
   <rotationref ref="rPMTRotation1"/>
  </physvol>
EOF
    }
  }
	print CRYOSTAT <<EOF;
 </volume>
</structure>
</gdml>
EOF

   close(CRYOSTAT);
}

#Generates Tia Miceli's scintillator veto wall
sub gen_vetoWall()
{
  # Set up the output file.
  $VW = "micro-vetoWall" . $suffix . ".gdml";
  push (@gdmlFiles, $VW); # Add file to list of GDML fragments
  $VW = ">" . $VW;
  open(VW) or die("Could not open file $VW for writing");
  
  print VW <<EOF;
  <?xml version='1.0'?>
  <gdml>
    <solids>
      <box lunit="cm" name="AuxDet0"    x="10"    y="20"    z="2" />
      <box lunit="cm" name="AuxDet1"    x="10"    y="20"    z="1" />
  
  </solids>
  <structure>
    <volume name="volAuxDet0">
      <materialref ref="Polystyrene"/>
      <solidref ref="AuxDet0"/>
    </volume>
    <volume name="volAuxDet1">
      <materialref ref="Polystyrene"/>
      <solidref ref="AuxDet1"/>
    </volume>
  </structure>
EOF
}

# Parameterize the cryostat's surroundings.
sub gen_enclosure()
{
    # Set up the output file.
    $GDML = "micro-enclosure" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="DetEnclosure" lunit="cm"
   x="$DetEnclosureWidth" y="$DetEnclosureHeight" z="$DetEnclosureLength"
 />
</solids>

<structure>
 <volume name="volDetEnclosure">
  <materialref ref="Air"/>
  <solidref ref="DetEnclosure"/>
  <physvol>
   <volumeref ref="volCryostat"/>
   <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
  </physvol>
EOF
#if the extra pieces in the enclosure are desired, put them in the correct places within volEnclosure
  if ( $enclosureExtras eq "on" ) {
    print GDML <<EOF;
      <physvol>
        <volumeref ref="volInsulation"/>
        <position name="posInsulation" unit="cm" x="0" y="0" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volPlatform"/>
        <position name="posPlatform" unit="cm" x="0" y="292.74" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volColumn"/>
        <position name="posColumn1" unit="cm" x="266" y="-121.261" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volColumn"/>
        <position name="posColumn2" unit="cm" x="-266" y="-121.261" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volTankBox1"/>
        <position name="posTank1_1" unit="cm" x="50" y="419" z="-600"/>
      </physvol>
      <physvol>
       <volumeref ref="volTankBox1"/>
       <position name="posTank1_2" unit="cm" x="-50" y="419" z="-600"/>
      </physvol>
      <physvol>
       <volumeref ref="volStandSubtraction"/>
       <position name="posStand1" unit="cm" x="0" y="-236" z="335.28"/>
      </physvol>
      <physvol>
         <volumeref ref="volStandSubtraction"/>
       <position name="posStand2" unit="cm" x="0" y="-236" z="-335.28"/>
      </physvol>
      <physvol>
       <volumeref ref="volStandConcrete"/>
       <position name="posConcreteStand1" unit="cm" x="0" y="-424" z="335.28"/>
      </physvol>
      <physvol>
         <volumeref ref="volStandConcrete"/>
       <position name="posConcreteStand2" unit="cm" x="0" y="-424" z="-335.28"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack1" unit="cm" x="-75" y="408.91" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack2" unit="cm" x="-75" y="408.91" z="93.85"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack3" unit="cm" x="-75" y="408.91" z="187.7"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack4" unit="cm" x="-75" y="408.91" z="275.15"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack5" unit="cm" x="-75" y="408.91" z="-142.15"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack6" unit="cm" x="-75" y="408.91" z="-243.8"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack7" unit="cm" x="-75" y="408.91" z="-326.9"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack8" unit="cm" x="127.9" y="408.91" z="-326.9"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack9" unit="cm" x="127.9" y="408.91" z="-387.7"/>
      </physvol>
      <physvol>
        <volumeref ref="volRack"/>
        <position name="posRack10" unit="cm" x="53.07" y="408.91" z="335.85"/>
        <rotationref ref="rPlus90AboutY"/>
      </physvol>
      <physvol>
         <volumeref ref="volFloorTankBox1"/>
         <position name="posfloortankbox1" unit="cm" x="-450" y="100+.001-530" z="335.28"/>
      </physvol>
      <physvol>
         <volumeref ref="volFloorTankBox1"/>
         <position name="posfloortankbox1_2" unit="cm" x="100" y="100+.001-530" z="-445"/>
      </physvol>
      <physvol>
         <volumeref ref="volFloorTankBox1"/>
         <position name="posfloortankbox1_3" unit="cm" x="100" y="100+.001-530" z="-521.5"/>
      </physvol>
      <physvol>
         <volumeref ref="volFloorTankBox1"/>
         <position name="posfloortankbox1_4" unit="cm" x="360" y="100+.001-530" z="-335.28"/>
      </physvol>
      <physvol>
         <volumeref ref="volFloorTankBox1"/>
         <position name="posfloortankbox1_5" unit="cm" x="360+60" y="100+.001-530" z="-335.28-60"/>
      </physvol>
      <physvol>
        <volumeref ref="volWalkway"/>
        <position name="posExtraPlatform" unit="cm" x="0" y="268" z="0"/>
      </physvol>
      <physvol>
        <volumeref ref="volWalkway"/>
        <position name="posWalkway" unit="cm" x="-212" y="181.74" z="0"/>
      </physvol>
      <physvol>
         <volumeref ref="volFloorTankBox2"/>
         <position name="posfloortankbox2_1" unit="cm" x="-20" y="100+.001-530" z="-500"/>
      </physvol>
      <physvol>
         <volumeref ref="volFloorTankBox2"/>
         <position name="posfloortankbox2_2" unit="cm" x="-136" y="40+80-530" z="-500"/>
      </physvol>
      <physvol>
         <volumeref ref="volPump"/>
         <position name="posPump" unit="cm" x="0" y="-424" z="200"/>
      </physvol>
EOF
}

  if ( $granite_block eq "on" ) {
    print GDML <<EOF;
     <!-- 06-Dec-2011 WGS: The x-position of the center of the granite block is =
              - CryostatRadius - GraniteBlockX/2 - room for scintillator modules
           The z-position is approximate; it will change when we know the
           position of other equipment in the hall. -->
     <physvol>
       <volumeref ref="volGraniteBlock"/>
       <position name="posGraniteBlock" unit="cm" x="-(193+(160/2)+10)" y="0" z="233"/>
     </physvol>
     <!-- 06-Dec-2011 WGS: The first two scintillator modules form an "X" on the
          low-x side of the granite block; the second two form an "X" on the far side. 
          I have no idea how these modules will be mounted; even though they're 
          2.52cm thick, I'm allowing 5cm for each module. -->
     <physvol>
       <volumeref ref="volDoubleChoozScintModule_0"/>
       <position name="posDoubleChoozScintModule_0" unit="cm" x="-(193+(5*0.5))" y="0" z="233" />
       <rotation name="rotDoubleChoozScintModule_0" unit="deg" x="45" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volDoubleChoozScintModule_1"/>
       <position name="posDoubleChoozScintModule_1" unit="cm" x="-(193+(5*1.5))" y="0" z="233" />
       <rotation name="rotDoubleChoozScintModule_1" unit="deg" x="135" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volDoubleChoozScintModule_2"/>
       <position name="posDoubleChoozScintModule_2" unit="cm" x="-(193+160+(5*2.5))" y="0" z="233" />
       <rotation name="rotDoubleChoozScintModule_2" unit="deg" x="45" y="0" z="0" />
     </physvol>
     <physvol>
       <volumeref ref="volDoubleChoozScintModule_3"/>
       <position name="posDoubleChoozScintModule_3" unit="cm" x="-(193+160+(5*3.5))" y="0" z="233" />
       <rotation name="rotDoubleChoozScintModule_3" unit="deg" x="135" y="0" z="0" />
     </physvol>
 
EOF
  }
  
  if ( $vetoWall_switch eq "on" ) {
    print GDML <<EOF;
    
    <physvol>
    <volumeref ref="volAuxDet0"/>
    <position name="posAuxDet0" unit="cm" x="0" y="0" z="-800"/>
    <rotation name="rAuxDet0" unit="deg" x="0" y="0" z="0"/>
    </physvol>
    <physvol>
    <volumeref ref="volAuxDet1"/>
    <position name="posAuxDet1" unit="cm" x="0" y="0" z="-805"/>
    <rotation name="rAuxDet1" unit="deg" x="0" y="0" z="0"/>
    </physvol>
    
EOF
  }

    print GDML <<EOF;
 </volume>
</structure>
</gdml>
EOF

   close(GDML);
}


# Parameterize the dirt mound that surrounds the enclosure.
sub gen_world()
{
    # Set up the output file.
    $GDML = "micro-world" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
  <box name="World" 
    lunit="cm" 
    x="$WorldWidth" 
    y="$WorldHeight" 
    z="$WorldLength"/>
  <tube name="Ground"
    rmin="310*2.54"
    rmax="((50*12)+310)*2.54"
    z="41*12*2.54"
    deltaphi="360" 
    lunit="cm"
    aunit="deg"/>
  <tube name="ConcreteEnclosure"
    rmin="292*2.54"
    rmax="310*2.54"
    z="38*12*2.54"
    deltaphi="360" 
    lunit="cm"
    aunit="deg"/>
  <tube name="PolystyreneEnclosure"
    rmin="310*2.54"
    rmax="312*2.54"
    z="(38*12+36)*2.54"
    deltaphi="360"
    lunit="cm"
    aunit="deg"/>
  <tube name="ConcreteEnclosureBottom"
    rmin="0"
    rmax="310*2.54"
    z="36*2.54"
    deltaphi="360" 
    lunit="cm"
    aunit="deg"/>
  <tube name="PolystyreneEnclosureBottom"
   rmax="292*2.54"
   z="2*2.54"
   deltaphi="360"
   lunit="cm"
   aunit="deg"/>
  <tube name="Overburden"
    rmin="0"
    rmax="584*2.54"
    z="10*12*2.54"
    deltaphi="360" 
    lunit="cm"
    aunit="deg"/>
</solids>

<structure>
  <volume name="volGround" >
    <materialref ref="Dirt" />
    <solidref ref="Ground" />
  </volume>
  <volume name="volOverburden" >
    <materialref ref="Dirt" />
    <solidref ref="Overburden" />
  </volume>
  <volume name="volPolystyreneEnclosure" >
    <materialref ref="Polystyrene" />
    <solidref ref="PolystyreneEnclosure" />
  </volume>
  <volume name="volConcreteEnclosure" >
    <materialref ref="Concrete" />
    <solidref ref="ConcreteEnclosure" />
  </volume>
  <volume name="volPolystyreneEnclosureBottom">
     <materialref ref="Polystyrene" />
     <solidref ref="PolystyreneEnclosureBottom"/>
  </volume>
  <volume name="volConcreteEnclosureBottom" >
    <materialref ref="Concrete" />
    <solidref ref="ConcreteEnclosureBottom" />
  </volume>
  <volume name="volWorld" >
    <materialref ref="Air"/> 
    <solidref ref="World"/>
    <physvol>
      <volumeref ref="volConcreteEnclosure"/>
      <position name="posConcreteEnclosure" unit="cm" x="0.5*$TPCActiveDepth" y="36*2.54/2" z="0.5*$TPCWirePlaneLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
    <physvol>
      <volumeref ref="volPolystyreneEnclosure"/>
      <position name="posPolystyreneEnclosure" unit="cm" x="0.5*$TPCActiveDepth" y="0" z="0.5*$TPCWirePlaneLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>  
    <physvol>
      <volumeref ref="volConcreteEnclosureBottom"/>
      <position name="posConcreteEnclosureBottom" unit="cm" x="0.5*$TPCActiveDepth" y="-38*12*2.54/2" z="0.5*$TPCWirePlaneLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>
   <physvol>
      <volumeref ref="volPolystyreneEnclosureBottom"/>
      <position name="posPolystyreneEnclosureBottom" unit="cm" x="0.5*$TPCActiveDepth" y="-(38*12 - 36)*2.54/2" z="0.5*$TPCWirePlaneLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol> 
    <physvol>
       <volumeref ref="volGround"/>
      <position name="posGround" unit="cm" x="0.5*$TPCActiveDepth" y="0" z="0.5*$TPCWirePlaneLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol>  
    <!--physvol>
      <volumeref ref="volOverburden"/>
      <position name="posOverburden" unit="cm" x="0.5*$TPCActiveDepth" y="(41-10)*12*2.54/2" z="0.5*$TPCWirePlaneLength"/>
      <rotationref ref="rPlus90AboutX"/>
    </physvol-->
    <physvol>
      <volumeref ref="volDetEnclosure"/>
      <position name="posDetEnclosure" unit="cm" x="0.5*$TPCActiveDepth" y="0" z="0.5*$TPCWirePlaneLength"/>
    </physvol>
  </volume> 
</structure>
</gdml>
EOF

   close(GDML);
}

sub gen_testbox()
{
    $test_switch="on";
    $GDML = "micro-testbox" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Define the solids and structures: the wires and the TPC wire plane.
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="TEST"
  lunit="cm"
  x="20"
  y="20"
  z="20"/>
</solids>
<structure>
 <volume name="volTEST">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="TEST"/>
 </volume>
</structure>
</gdml>
EOF

    close(GDML);
}
####################################
sub gen_enclosureExtras()
{
    $GDML = "micro-enclosureExtras" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Define the solids and structures associated with the "enclosureExtras" (platform, two walkways, insulation, cement stands and foam saddles, tanks on the floor  (7 - 5 big, 2, small) and platform (2), computer racks on top of the platform (10)). Names of particular volumes are explanatory
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
<tube name="Insulation"
  rmin="$CryostatOuterRadius+0.1"
  rmax="$CryostatOuterRadius+0.1+$inch*$InsulationThickness"
  z="$CryostatLength"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="Platform"
  x="550"
  y="18"
  z="1500"
  lunit="cm"/>
<box name="Column"
  x="16.79"
  y="798"
  z="16.79"
  lunit="cm"/>
<tube name="Tank1"
  rmin="30.35"
  rmax="38"  
  z="140.3" 
  deltaphi="360" 
  aunit="deg" 
  lunit="cm"/> 
<tube name="TankCap1"
  rmax="38"
  z="40"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<box name="TankBox1"
  x="76.001"
  y="221"
  z="76.001"
  lunit="cm"/>

<box name="standBox"
  x="507"
  y="196.27"
  z="91.44"
  lunit="cm"/>
<tube name="standSubTube"
  rmax="233.75"
  z="125"
  deltaphi="360"
  aunit="deg"
  lunit="cm"/>
<subtraction name="standSubtraction">
  <first ref="standBox"/>
  <second ref="standSubTube"/>
  <position name="posStandSubtract" unit="cm" x="0" y="235.92" z="0"/>
</subtraction>
<box name="standConcrete"
  x="548.64"
  y="182.89"
  z="91.44"
  lunit="cm"/>

<box name="rackBox"
  x="$RackX"
  y="$RackY"
  z="$RackZ"
  lunit="cm"/>
<box name="rackX"
  x="$RackX"
  y="$RackThickness"
  z="$RackThickness"
  lunit="cm"/>
<box name="rackY"
  x="$RackThickness"
  y="$RackY-2*$RackThickness-0.001"
  z="$RackThickness"
  lunit="cm"/>
<box name="rackZ"
  x="$RackThickness"
  y="$RackThickness"
  z="$RackZ-2*$RackThickness-0.001"
  lunit="cm"/>

<box name="floorTankBox1"
   x="59.7"
   y="200"
   z="59.7"
   lunit="cm"/>
<tube name="floorTank1"
   rmin="45.4/2"
   rmax="59.6/2"
   z="159.8"
   deltaphi="360"
   aunit="deg"
   lunit="cm"/>
<tube name="floorTankCap1"
   rmax="59.6/2"
   z="20.0"
   deltaphi="360"
   aunit="deg"
   lunit="cm"/>
<box name="floorTankBox2"
   x="52.2"
   y="67"
   z="52.2"
   lunit="cm"/>
<tube name="floorTank2"
   rmin="40.0/2"
   rmax="52.1/2"
   z="53.0"
   deltaphi="360"
   aunit="deg"
   lunit="cm"/>
<tube name="floorTankCap2"
   rmax="52.1/2"
   z="52.1/4"
   deltaphi="360"
   aunit="deg"
   lunit="cm"/>

<box name="Walkway"
   x="91"
   y="24"
   z="1200"
   lunit="cm"/>

<box name="Pump"
   x="100"
   y="196"
   z="91.7"
   lunit="cm"/>
</solids>
<structure>
  <volume name="volInsulation">
     <materialref ref="PU_foam_light"/>
     <solidref ref="Insulation"/>
  </volume>
    <volume name="volPlatform">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="Platform"/>
    </volume>
    <volume name="volColumn">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="Column"/>
    </volume>
    <volume name="volTankCap1">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="TankCap1"/>
    </volume>
    <volume name="volTank1">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="Tank1"/>
    </volume>
    <volume name="volStandSubtraction">
      <materialref ref="PU_foam_dense"/>
      <solidref ref="standSubtraction"/>
    </volume>

    <volume name="volTankBox1">
      <materialref ref="Air"/>
      <solidref ref="TankBox1"/>
      <physvol>
        <volumeref ref="volTank1"/>
        <position name="posvolTank1" unit="cm" x="0" y="0" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volTankCap1"/> 
        <position name="posvolTankCap1_1" unit="cm" x="0" y="90.15+.001" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volTankCap1"/> 
        <position name="posvolTankCap1_2" unit="cm" x="0" y="-90.15-.001" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
    </volume>

    <volume name="volStandConcrete">
      <materialref ref="Concrete"/>
      <solidref ref="standConcrete"/>
    </volume>

    <volume name="volRackX">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="rackX"/>
    </volume>
    <volume name="volRackY">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="rackY"/>
    </volume>
    <volume name="volRackZ">
      <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
      <solidref ref="rackZ"/>
    </volume>

    <volume name="volRack">
      <materialref ref="Air"/>
      <solidref ref="rackBox"/>
      <physvol>
        <volumeref ref="volRackX"/>
        <position name="posRackX1" unit="cm" x="0" y="($RackY-$RackThickness)/2" z="($RackZ-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackX"/>
        <position name="posRackX2" unit="cm" x="0" y="-($RackY-$RackThickness)/2" z="($RackZ-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackX"/>
        <position name="posRackX3" unit="cm" x="0" y="($RackY-$RackThickness)/2" z="-($RackZ-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackX"/>
        <position name="posRackX4" unit="cm" x="0" y="-($RackY-$RackThickness)/2" z="-($RackZ-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackY"/>
        <position name="posRackY1" unit="cm" y="0" x="($RackX-$RackThickness)/2" z="($RackZ-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackY"/>
        <position name="posRackY2" unit="cm" y="0" x="-($RackX-$RackThickness)/2" z="($RackZ-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackY"/>
        <position name="posRackY3" unit="cm" y="0" x="($RackX-$RackThickness)/2" z="-($RackZ-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackY"/>
        <position name="posRackY4" unit="cm" y="0" x="-($RackX-$RackThickness)/2" z="-($RackZ-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackZ"/>
        <position name="posRackZ1" unit="cm" z="0" y="($RackY-$RackThickness)/2" x="($RackX-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackZ"/>
        <position name="posRackZ2" unit="cm" z="0" y="($RackY-$RackThickness)/2" x="-($RackX-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackZ"/>
        <position name="posRackZ3" unit="cm" z="0" y="-($RackY-$RackThickness)/2" x="($RackX-$RackThickness)/2"/>
      </physvol>
      <physvol>
        <volumeref ref="volRackZ"/>
        <position name="posRackZ4" unit="cm" z="0" y="-($RackY-$RackThickness)/2" x="-($RackX-$RackThickness)/2"/>
      </physvol>  
    </volume>

  <volume name="volFloorTank1">
	<materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
	<solidref ref="floorTank1"/>
  </volume>

  <volume name="volFloorTankCap1">
	<materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
	<solidref ref="floorTankCap1"/>
  </volume>

  <volume name="volFloorTank2">
	<materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
	<solidref ref="floorTank2"/>
  </volume>

  <volume name="volFloorTankCap2">
	<materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
	<solidref ref="floorTankCap2"/>
  </volume>

 <volume name="volFloorTankBox1">
      <materialref ref="Air"/>
      <solidref ref="floorTankBox1"/>
      <physvol>
        <volumeref ref="volFloorTank1"/>
        <position name="posFloorTank1" unit="cm" x="0" y="0" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volFloorTankCap1"/> 
        <position name="posFloorTankCap1_1" unit="cm" x="0" y="89.9001" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volFloorTankCap1"/> 
        <position name="posFloorTankCap1_2" unit="cm" x="0" y="-89.9001" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
    </volume>

    <volume name="volFloorTankBox2">
      <materialref ref="Air"/>
      <solidref ref="floorTankBox2"/>
      <physvol>
        <volumeref ref="volFloorTank2"/>
        <position name="posFloorTank2" unit="cm" x="0" y="6.033" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
      <physvol>
        <volumeref ref="volFloorTankCap2"/>
        <position name="posFloorTankCap2" unit="cm" x="0" y="-26.98" z="0"/>
        <rotationref ref="rPlus90AboutX"/>
      </physvol>
    </volume>

     <volume name="volWalkway">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="Walkway"/>
     </volume>


     <volume name="volPump">
        <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
        <solidref ref="Pump"/>
     </volume>

</structure>
</gdml>
EOF

    close(GDML);
}
#####################################3


sub write_fragments()
{
    # The output file is a list of the GDML sub-files created by this
    # script.

    if ( ! defined $output )
    {
	$output = "-"; # write to STDOUT 
    }

    # Set up the output file.
    $OUTPUT = ">" . $output;
    open(OUTPUT) or die("Could not open file $OUTPUT");

    print OUTPUT <<EOF;
<?xml version='1.0'?>

<!-- Input to Geometry/gdml/make_gdml.pl; define the GDML fragments
     that will be zipped together to create a detector description. 
     -->

<config>

   <constantfiles>

      <!-- These files contain GDML <constant></constant>
           blocks. They are read in separately, so they can be
           interpreted into the remaining GDML. See make_gdml.pl for
           more information. 
	   -->
	   
EOF

    foreach $filename (@defFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </constantfiles>

   <gdmlfiles>

      <!-- The GDML file fragments to be zipped together. -->

EOF

    foreach $filename (@gdmlFiles)
    {
	print OUTPUT <<EOF;
      <filename> $filename </filename>
EOF
    }

    print OUTPUT <<EOF;

   </gdmlfiles>

</config>
EOF

    close(OUTPUT);
}
