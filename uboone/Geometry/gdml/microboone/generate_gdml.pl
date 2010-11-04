#!/usr/bin/perl

# This program creates GDML sub-files, with values supplied by user
# parameters.  Geometry/gdml/make_gdml.pl "zips" together those
# sub-files to make a single detector description.

# Packages
use Math::Trig;
use XML::LibXML;
use Getopt::Long;

# Get the input parameters from an XML file. Optionally append a
# suffix to the GDML sub-files we create.

GetOptions( "input|i:s" => \$input,
	    "help|h" => \$help,
	    "suffix|s:s" => \$suffix,
	    "output|o:s" => \$output);

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
    # strings in them (like "kInch").

    eval "\$$name = '$value'";
}

# Our calculations and constants depend on the geometry of the wires.
$SinUVAngle = sin( deg2rad($UVAngle) );
$CosUVAngle = cos( deg2rad($UVAngle) );
$TanUVAngle = tan( deg2rad($UVAngle) );

# The routines that create the GDML sub-files. Most of the explanatory
# comments are in gen_defs().

gen_defs();
gen_rotations();
gen_materials();
# following two subroutines generate wireplanes which slow down visualization, so for the sake
# of developing this without waiting too long in between, i'm commenting these out
#if ( $NumberOfTPCPlanes == 3 ) { gen_microvertplane(); }
#gen_microplane();
gen_tpc();
gen_cryostat();
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
    # Set up the output file.
    $CONSTANTS = "micro-defs" . $suffix . ".gdml";
    push (@defFiles, $CONSTANTS); # Add file to list of constant files
    $CONSTANTS = ">" . $CONSTANTS;
    open(CONSTANTS) or die("Could not open file $CONSTANTS for writing");

    # Create some math constants.
    my $pi = pi;

    # Though it's not strictly necessary, make each sub-file a valid
    # XML (if not GDML) document; that makes it accessible to an XML
    # parser if needed.

    # Here is a neat way to print out a block of text without getting
    # involved in a lot of messy quoting and formatting with print
    # statements.

    print CONSTANTS <<EOF;
<?xml version='1.0'?>
<define>
<constant name="kInch"	value="2.54" />
<constant name="kPi"	value="$pi" />
<constant name="kDetEnclosureWidth"	  value="$DetEnclosureWidth" />
<constant name="kDetEnclosureHeight"	  value="$DetEnclosureHeight" />
<constant name="kDetEnclosureLength"	  value="$DetEnclosureLength" />
<constant name="kDirtThickness"           value="$DirtThickness" />
<constant name="kWorldW"                  value="100.0*kDetEnclosureWidth"/>
<constant name="kWorldH"                  value="100.0*kDetEnclosureHeight"/>
<constant name="kWorldL"                  value="100.0*kDetEnclosureLength"/>

<constant name="kTPCWidth"		  value="$TPCWidth" />
<constant name="kTPCLength"		  value="$TPCLength" />
<constant name="kTPCDepth"		  value="$TPCDepth" />
<constant name="kTPCWallThickness"        value="$TPCWallThickness" />

<constant name="kTPCWirePlaneThickness"   value="$TPCWirePlaneThickness" />
<constant name="kTPCWireThickness"	  value="$TPCWireThickness" />
<constant name="kTPCWirePlaneWidth"       value="kTPCWidth" />
<constant name="kTPCWirePlaneLength"	  value="kTPCLength" />
<constant name="kTPCWirePitch"            value="$TPCWirePitch"/>
<constant name="kSinUVAngle"              value="$SinUVAngle"/>
<constant name="kCosUVAngle"              value="$CosUVAngle"/>
<constant name="kTanUVAngle"              value="$TanUVAngle"/>
<constant name="kTPCWireXPitch"           value="kTPCWirePitch/kCosUVAngle"/>
<constant name="kCryoTPCMountLength" 	  value="40">
<constant name="kCryoTPCMountWidth" 	  value="20">
<constant name="kCryoTPCMountHeight" 	  value="20">
<constant name="kRailLength"               value="1200">
<constant name="kRailWidth"                value="4">
<constant name="kRailHeight"               value="8">
</define>
EOF

   close(CONSTANTS);
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
   <rotation name="rPlusUVAngleAboutX"	unit="deg" x="$WirePlusRotation" y="0"   z="0"/>
   <rotation name="rPlus180AboutX"	unit="deg" x="180" y="0"   z="0"/>
   <rotation name="rMinusUVAngleAboutX" unit="deg" x="$WireMinusRotation" y="0"   z="0"/>
   <rotation name="rMinus90AboutX" unit="deg" x="-90" y="0"   z="0"/>
   <rotation name="rPlus30AboutY"  unit="deg" x="0"   y="30"  z="0"/>
   <rotation name="rPlus60AboutY"  unit="deg" x="0"   y="60"  z="0"/>
   <rotation name="rPlus90AboutY"  unit="deg" x="0"   y="90"  z="0"/>
   <rotation name="rPlus180AboutY" unit="deg" x="0"   y="180" z="0"/>
   <rotation name="rMinus90AboutY" unit="deg" x="0"   y="-90" z="0"/>
   <rotation name="rPlus90AboutZ"  unit="deg" x="0"   y="0"   z="90"/>
   <rotation name="rMinus90AboutZ"  unit="deg" x="0"   y="0"   z="-90"/>
   <rotation name="rPlus180AboutZ"	unit="deg" x="0"   y="0"   z="180"/>
   <rotation name="rMinus180AboutZ"	unit="deg" x="0"   y="0"   z="-180"/>
   <rotation name="rPlus90AboutXPlus90AboutZ"  unit="deg" x="90" y="0"   z="90"/>
   <rotation name="rPlus90AboutXPlus180AboutZ" unit="deg" x="90" y="0"   z="180"/>
   <rotation name="rPlus90AboutYMinus90AboutZ" unit="deg" x="0"  y="90"  z="0"/>
   <rotation name="rPlus90AboutXMinus90AboutY" unit="deg" x="90" y="-90" z="0"/>
   <rotation name="rPlus90AboutXMinus90AboutZ" unit="deg" x="90" y="0"   z="-90"/>
   <rotation name="rPlus90AboutXPlus90AboutY"  unit="deg"  x="90" y="90" z="0"/>
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
    my $NumberWires = int( $TPCLength / $TPCWirePitch ) - 1;

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
  rmax="0.5*kTPCWireThickness"
  z="kTPCWidth"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
<box name="TPCPlaneVert"
  x="kTPCWirePlaneThickness" 
  y="kTPCWidth" 
  z="kTPCLength"
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
	     <position name="posTPCWireVert$i" unit="cm" z="-0.5*kTPCLength+kTPCWirePitch*($i+1)" x="0" y="0"/>
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
    my $NumberWires = $TPCLength / $TPCWirePitch - 1;

    $GDML = "micro-plane" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    # Calculate the number of wire ends on a given y-edge of the plane.
    my $TPCYWirePitch = $TPCWirePitch / $CosUVAngle;
    my $NumberWiresPerEdge = int( $TPCLength / $TPCYWirePitch );

    # How many side wires will be "cut off" by the lower or higher
    # z-edge?
    my $NumberSideWires = int( $TanUVAngle * $TPCWidth / $TPCYWirePitch );

    # The number of full-length "center" wires.
    my $NumberCenterWires = $NumberWiresPerEdge - $NumberSideWires;

    # define the solids
    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
EOF

    # wires on either end of the tpc
    for($i = 0; $i < $NumberSideWires; ++$i)
    {
	print GDML <<EOF;
<tube name="TPCWire$i"
  rmax="0.5*kTPCWireThickness"
  z="kTPCWireXPitch*($i+1)/kSinUVAngle"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
EOF
    }

    # The solids for the middle wire and the TPC wire plane, and start off the structures.
    print GDML <<EOF;
<tube name="TPCWireCommon"
  rmax="0.5*kTPCWireThickness"
  z="kTPCWidth/kCosUVAngle"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
<box name="TPCPlane"
  x="kTPCWirePlaneThickness"
  y="kTPCWidth"
  z="kTPCLength"
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
      print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWire$i"/> 
     <position name="posTPCWire$i" unit="cm" y="-0.5*kTPCWidth+0.5*($i+1)*kTPCWireXPitch/kTanUVAngle" z="-0.5*kTPCLength+0.5*kTPCWireXPitch*($i+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  }

  # The wires in the middle.
  for ($i = 0; $i < $NumberCenterWires; ++$i)
  {
      my $j = $NumberSideWires+$i;
      print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWireCommon"/>
     <position name="posTPCWire$j" unit="cm" y="0" z="-0.5*kTPCWirePlaneLength+kTPCWireXPitch*(0.5*$NumberSideWires + $i+1)" x="0"/>
     <rotationref ref="rPlusUVAngleAboutX"/>
    </physvol>
EOF
  }

  # the wires at the +z end
  for ($i = 0; $i < $NumberSideWires; ++$i)
  {
      my $j = $NumberSideWires-$i-1;
      my $k = $NumberCenterWires+$NumberSideWires+$i;

      print GDML <<EOF;
    <physvol>
     <volumeref ref="volTPCWire$j"/>
     <position name="posTPCWire$k" unit="cm" y="0.5*kTPCWidth-0.5*($j+1)*kTPCWireXPitch/kTanUVAngle" z="0.5*kTPCWirePlaneLength-0.5*kTPCWireXPitch*($j+1)" x="0"/>
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


# Parameterize the TPC and the planes within it.
sub gen_tpc()
{
    # Set up the output file.
    $GDML = "micro-tpc" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <box name="TPCBackWall"
  lunit="cm"
  z="kTPCLength"
  x="kTPCWallThickness"
  y="kTPCWidth"/>
 <box name="TPCVertWall"
  lunit="cm"
  x="kTPCDepth+kTPCWallThickness"
  y="kTPCWidth"
  z="kTPCWallThickness"/>
 <box name="TPCHorizWall"
  lunit="cm"
  x="kTPCDepth+kTPCWallThickness"
  y="kTPCWallThickness"
  z="kTPCLength+2*kTPCWallThickness"/>
 <box name="TPC"
  lunit="cm"
  x="kTPCDepth+kTPCWallThickness+2*kTPCWirePlaneThickness"
  y="kTPCWidth+2*kTPCWallThickness"
  z="kTPCLength+2*kTPCWallThickness"/>
</solids>
<structure>
 <volume name="volTPCBackWall">
  <materialref ref="G10"/>
  <solidref ref="TPCBackWall"/>
 </volume>
 <volume name="volTPCHorizWall">
  <materialref ref="G10"/>
  <solidref ref="TPCHorizWall"/>
 </volume>
 <volume name="volTPCVertWall">
  <materialref ref="G10"/>
  <solidref ref="TPCVertWall"/>
 </volume>
 <volume name="volTPC">
  <materialref ref="LAr"/>
  <solidref ref="TPC"/>
  <physvol>
   <volumeref ref="volTPCBackWall"/>
   <position name="posTPCBackWall" unit="cm" x="-0.5*kTPCDepth" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCVertWall"/>
   <position name="posTPCVertWallN" unit="cm" x="0" y="0" z="0.5*(kTPCLength+kTPCWallThickness)"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCVertWall"/>
   <position name="posTPCVertWallS" unit="cm"  x="0" y="0" z="-0.5*(kTPCLength+kTPCWallThickness)"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCHorizWall"/>
   <position name="posTPCBottomWall" unit="cm" x="0" y="-0.5*(kTPCWidth+kTPCWallThickness)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCHorizWall"/>
   <position name="posTPCTopWall" unit="cm" x="0" y="0.5*(kTPCWidth+kTPCWallThickness)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCPlane"/>
   <position name="posTPCPlane1" unit="cm" x="0.5*kTPCDepth" y="$UPlaneYOffset" z="$UPlaneZOffset"/>
   <rotationref ref="rPlus180AboutZ"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPCPlane"/>
   <position name="posTPCPlane2" unit="cm" x="0.5*kTPCDepth+1.0*kTPCWirePlaneThickness" y="$VPlaneYOffset" z="$VPlaneZOffset"/>
  </physvol>
EOF

   if ( $NumberOfTPCPlanes == 3 )
   {
       print GDML <<EOF;
  <physvol>
   <volumeref ref="volTPCPlaneVert"/>
   <position name="posTPCPlaneVert" unit="cm" x="0.5*kTPCDepth+2.0*kTPCWirePlaneThickness" y="0" z="$YPlaneZOffset"/>
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


# Parameterize the steel cryostat that encloses the TPC.
sub gen_cryostat()
{
    # Set up the output file.
    $GDML = "micro-cryostat" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    print GDML <<EOF;
<?xml version='1.0'?>
<gdml>
<solids>
 <tube name="Cryostat" 
  rmax="$CryostatOuterRadius"
  z="$CryostatLength+2*$CryostatEndcapThickness"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
<tube name="SteelTube"
  rmin="$CryostatInnerRadius"
  rmax="$CryostatOuterRadius"
  z="$CryostatLength"
  deltaphi="2*kPi"
  aunit="rad"
  lunit="cm"/>
<box name="CryoTPCRail"
  x="10" 
  y="10" 
  z="1200"
  lunit="cm"/>
<box name="CryoPlaneMount"
  x="10"
  y="10"
  z="10"
  lunit="cm"/>
</solids>

<structure>
 <volume name="volSteelTube">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="SteelTube"/>
 </volume>
 <volume name="volCryoRailLeft">
  <materialref ref="ALUMINUM_Al"/>
  <solidref ref="CryoTPCRail"/>
 </volume>
 <volume name="volCryoRailRight">
  <materialref ref="ALUMINUM_Al"/>
  <solidref ref="CryoTPCRail"/>
 </volume>
 <volume name="volPlaneMountCathodeFront">
  <materialref ref="ALUMINUM_Al"/>
  <solidref ref="CryoPlaneMount"/>
 </volume>
 <volume name="volPlaneMountCathodeBack">
  <materialref ref="ALUMINUM_Al"/>
  <solidref ref="CryoPlaneMount"/>
 </volume>
 <volume name="volPlaneMountWireFront">
  <materialref ref="ALUMINUM_Al"/>
  <solidref ref="CryoPlaneMount"/>
 </volume>
 <volume name="volPlaneMountWireBack">
  <materialref ref="ALUMINUM_Al"/>
  <solidref ref="CryoPlaneMount"/>
 </volume>
 <volume name="volCryostat">
  <materialref ref="LAr"/>
  <solidref ref="Cryostat"/>
  <physvol>
   <volumeref ref="volSteelTube"/>
   <position name="posSteelTube" unit="cm" x="0" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volTPC"/>
   <position name="posTPC" x="0.0" y="0" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volCryoRailLeft"/>
   <position name="posCryoRailLeft" unit="cm" x="-30" y="-235.557321" z="100"/>
  </physvol>
  <physvol>
   <volumeref ref="volCryoRailRight"/>
   <position name="posCryoRailRight" unit="cm" x="30" y="-235.557321" z="100"/>
  </physvol>
  <physvol>
   <volumeref ref="volPlaneMountCathodeFront"/>
   <position name="posTPC" x="-200" y="-0.5*kTPCWirePlaneWidth" z="-0.25*kTPCWirePlaneLength"/>
  </physvol>
  <physvol>
   <volumeref ref="volPlaneMountCathodeBack"/>
   <position name="posTPC" x="-200" y="-0.5*kTPCWirePlaneWidth" z="0.25*kTPCWirePlaneLength"/>
  </physvol>
  <physvol>
   <volumeref ref="volPlaneMountWireFront"/>
   <position name="posTPC" x="200" y="-0.5*kTPCWirePlaneWidth" z="-0.25*kTPCWirePlaneLength"/>
  </physvol>
  <physvol>
   <volumeref ref="volPlaneMountWireBack"/>
   <position name="posTPC" x="200" y="-0.5*kTPCWirePlaneWidth" z="0.25*kTPCWirePlaneLength"/>
  </physvol>
 </volume>
</structure>
</gdml>
EOF

   close(GDML);
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
   x="kDetEnclosureWidth" y="kDetEnclosureHeight" z="kDetEnclosureLength"
 />
</solids>

<structure>
 <volume name="volDetEnclosure">
  <materialref ref="Air"/>
  <solidref ref="DetEnclosure"/>
  <physvol>
   <volumeref ref="volCryostat"/>
   <position name="posCryostat" unit="cm" x="0" y="0" z="0"/>
   <rotationref ref="rMinus180AboutZ"/>
  </physvol>
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
   <box name="World" lunit="cm" x="kWorldW" y="kWorldH" z="kWorldL"/>
   <box name="MoundTop" x="kDetEnclosureWidth" y="kDirtThickness" z="kDetEnclosureLength" lunit="cm" />
   <xtru name="MoundEW" lunit="cm" >
     <twoDimVertex x="0"                y="-0.5*kDetEnclosureHeight" />
     <twoDimVertex x="4*kDirtThickness" y="-0.5*kDetEnclosureHeight" />
     <twoDimVertex x="kDirtThickness"   y="0.5*kDetEnclosureHeight+kDirtThickness" />
     <twoDimVertex x="0" y="0.5*kDetEnclosureHeight+kDirtThickness" />
     <section zOrder="0" zPosition="-0.5*kDetEnclosureLength" xOffset="0" yOffset="0" scalingFactor="1" />
     <section zOrder="1" zPosition="0.5*kDetEnclosureLength"  xOffset="0" yOffset="0" scalingFactor="1" />  
   </xtru>
   <xtru name="MoundNS" lunit="cm" >
     <twoDimVertex x="0"                y="-0.5*kDetEnclosureHeight" />
     <twoDimVertex x="4*kDirtThickness" y="-0.5*kDetEnclosureHeight" />
     <twoDimVertex x="kDirtThickness"   y="0.5*kDetEnclosureHeight+kDirtThickness" />
     <twoDimVertex x="0" y="0.5*kDetEnclosureHeight+kDirtThickness" />
     <section zOrder="0" zPosition="-0.5*kDetEnclosureWidth" xOffset="0" yOffset="0" scalingFactor="1" />
     <section zOrder="1" zPosition="0.5*kDetEnclosureWidth"  xOffset="0" yOffset="0" scalingFactor="1" />
   </xtru>
   <cone name="MoundCorner" lunit="cm" aunit="deg"
     rmin1="0" rmin2="0" rmax1="4*kDirtThickness" rmax2="kDirtThickness"
     z="kDetEnclosureHeight+kDirtThickness"
     startphi="90" deltaphi="90"
   />
</solids>

<structure>
 <volume name="volMoundTop" >
  <materialref ref="Dirt" />
  <solidref ref="MoundTop" />
 </volume>
 <volume name="volMoundEast" >
  <materialref ref="Dirt" />
  <solidref ref="MoundEW" />
 </volume>
 <volume name="volMoundWest" >
  <materialref ref="Dirt" />
  <solidref ref="MoundEW" />
 </volume>
 <volume name="volMoundNorth" >
  <materialref ref="Dirt" />
  <solidref ref="MoundNS" />
 </volume>
 <volume name="volMoundCorner" >
  <materialref ref="Dirt" />
  <solidref ref="MoundCorner" />
 </volume>
 <volume name="volWorld" >
  <materialref ref="Air"/> <solidref ref="World"/>
   <physvol>
    <volumeref ref="volMoundTop"/>
    <position name="posMoundTop" unit="cm" x="0.5*kTPCDepth" y="0.5*(kDetEnclosureHeight+kDirtThickness)" z="0.5*kTPCLength"/>
   </physvol>
   <physvol>
    <volumeref ref="volMoundCorner"/>
    <position name="posMoundCornerSW" unit="cm" x="-0.5*(kDetEnclosureWidth-kTPCDepth)" y="0.5*kDirtThickness" z="-0.5*(kDetEnclosureLength-kTPCLength)"/>
    <rotationref ref="rPlus90AboutX"/>
   </physvol>
   <physvol>
    <volumeref ref="volMoundCorner"/>
    <position name="posMoundCornerNW" unit="cm" x="0.5*(kDetEnclosureWidth+kTPCDepth)" y="0.5*kDirtThickness" z="0.5*(kDetEnclosureLength+kTPCLength)"/>
    <rotationref ref="rPlus90AboutXPlus180AboutZ"/>
   </physvol>
   <physvol>
    <volumeref ref="volMoundCorner"/>
    <position name="posMoundCornerNE" unit="cm" x="-0.5*(kDetEnclosureWidth-kTPCDepth)" y="0.5*kDirtThickness" z="0.5*(kDetEnclosureLength+kTPCLength)"/>
    <rotationref ref="rPlus90AboutXMinus90AboutZ"/>
   </physvol>
   <physvol>
    <volumeref ref="volMoundCorner"/>
    <position name="posMoundCornerSE" unit="cm" x="0.5*(kDetEnclosureWidth+kTPCDepth)" y="0.5*kDirtThickness" z="-0.5*(kDetEnclosureLength-kTPCLength)"/>
    <rotationref ref="rPlus90AboutXPlus90AboutZ"/>
   </physvol>
   <physvol>
    <volumeref ref="volMoundEast"/>
    <position name="posMoundEast" unit="cm" x="0.5*(kTPCDepth+kDetEnclosureWidth)" y="0" z="0.5*kTPCLength"/>
   </physvol>
   <physvol>
    <volumeref ref="volMoundWest"/>
    <position name="posMoundWest" unit="cm" x="-0.5*(kDetEnclosureWidth-kTPCDepth)" y="0" z="0.5*kTPCLength"/>
    <rotationref ref="rPlus180AboutY"/>
   </physvol>
   <physvol>
    <volumeref ref="volMoundNorth"/>
    <position name="posMoundNorth" unit="cm" x="0.5*kTPCDepth" y="0" z="0.5*(kTPCLength+kDetEnclosureLength)"/>
    <rotationref ref="rPlus90AboutY"/>
   </physvol>
   <physvol>
    <volumeref ref="volDetEnclosure"/>
    <position name="posDetEnclosure" unit="cm" x="0.5*kTPCDepth" y="0" z="0.5*kTPCLength"/>
   </physvol>
 </volume> 
</structure></gdml>
EOF

   close(GDML);
}

sub gen_fieldmesh() {
#this subroutine will produce general gdml fragment for the field mesh

    # Set up the output file.
    $GDML = "micro-fieldmesh" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    #put the GDML IN B/W EOF'S
    print GDML <<EOF;

EOF

   close(GDML);
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
<define>
 <constant name="kGroundPlateWidth"  	value="224" /> 
 <constant name="kGroundPlateHeight"    	value="0.1" /> 
 <constant name="kGroundPlateLength"       	value="1100" /> 
 <constant name="kGroundBeamWidth"  	value="2.5" /> 
 <constant name="kGroundBeamHeight"    	value="2.5" /> 
 <constant name="kGroundBeamLength"       	value="kGroundPlateLength" />
 <constant name="kGroundBeamThickness"  	value=".15" /> 
</define>
<solids>
 <box name="GroundPlane"
  lunit="cm"
  x="kGroundPlateWidth+0.5"
  y="kGroundPlateHeight+kGroundBeamHeight"
  z="kGroundPlateLength"/>
 <box name="GroundPlate"
  lunit="cm"
  x="kGroundPlateWidth"
  y="kGroundPlateHeight"
  z="kGroundPlateLength"/>
 <box name="GroundBeam"
  lunit="cm"
  x="kGroundBeamWidth"
  y="kGroundBeamHeight"
  z="kGroundBeamLength"/>
 <box name="GroundBeamIn"
  lunit="cm"
  x="kGroundBeamWidth-kGroundBeamThickness"
  y="kGroundBeamHeight-kGroundBeamThickness"
  z="kGroundBeamLength"/>
</solids>
<structure>
 <volume name="volGroundPlate">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="GroundPlate"/>
 </volume>
 <volume name="volGroundBeam">
  <materialref ref="STEEL_STAINLESS_Fe7Cr2Ni"/>
  <solidref ref="GroundBeam"/>
 </volume>
 <volume name="volGroundPlane">
  <materialref ref="LAr"/>
  <solidref ref="GroundPlane"/>
  <physvol>I
   <volumeref ref="volGroundPlate"/>
   <position name="posGroundPlate" unit="cm" x="0.5" y="-0.5*(kGroundBeamHeight)" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volGroundBeam"/>
   <position name="posGroundBeam1" unit="cm" x="-(0.5*(kGroundPlateWidth-kGroundBeamWidth)+0.5)" y="0.5*kGroundPlateHeight" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volGroundBeam"/>
   <position name="posGroundBeam2" unit="cm" x="-12" y="0.5*kGroundPlateHeight" z="0"/>
  </physvol>
  <physvol>
   <volumeref ref="volGroundBeam"/>
   <position name="posGroundBeam3" unit="cm" x="0.5*(kGroundPlateWidth-kGroundBeamWidth)-9" y="0.5*kGroundPlateHeight" z="0"/>
  </physvol>
 </volume> 
</structure>
</gdml>
EOF

    close(GROUNDPLATE);
}

sub gen_fieldcage() {
#this subroutine will produce general gdml fragment for the field cage
    # Set up the output file.
    $GDML = "micro-fieldcage" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    #put the GDML IN B/W EOF'S
    print GDML <<EOF;

EOF


   close(GDML);
}

sub gen_pmtrack() {
#this subroutine will produce general gdml fragment for the pmt rack setup
    # Set up the output file.
    $GDML = "micro-pmtrack" . $suffix . ".gdml";
    push (@gdmlFiles, $GDML); # Add file to list of GDML fragments
    $GDML = ">" . $GDML;
    open(GDML) or die("Could not open file $GDML for writing");

    #put the GDML IN B/W EOF'S
    print GDML <<EOF;

EOF

    close(GDML);
}


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

