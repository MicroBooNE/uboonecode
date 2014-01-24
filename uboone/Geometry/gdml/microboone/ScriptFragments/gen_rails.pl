#!/usr/bin/perl

$CryostatInnerRadius=190.46;
$CryostatLength=1100;
$x_pos=25;
$x_i=$x_pos;
$CryostatRailOffset=150;

$CryoRailLength=$CryostatLength-$CryostatRailOffset;

$CRYOSTAT = "micro-cryostat2" . $suffix . ".gdml";
push (@gdmlFiles, $CRYOSTAT); # Add file to list of CRYOSTAT fragments
$CRYOSTAT = ">" . $CRYOSTAT;
open(CRYOSTAT) or die("Could not open file $CRYOSTAT for writing");

for ( $j=-1; $j<2; $j+=2 ) 
{
	$x_pos=$x_i;
	$x_pos*=$j;
	print CRYOSTAT <<EOF;
 <xtru name="CryoRail$j" lunit="cm" >	
EOF

	for ( $i=0; $i<7; $i++ )
	{
		$y_pos=sqrt( $CryostatInnerRadius**2 - $x_pos**2 );
		print CRYOSTAT <<EOF;
  <twoDimVertex x="$x_pos" y="$y_pos" />
EOF
		$x_pos+=$j;
	}

	$x_f=$x_i*$j;
	print CRYOSTAT <<EOF;
  <twoDimVertex x="$x_f" y="$y_pos" />
  <section zOrder="0" zPosition="($CryostatLength/2)-$CryostatRailOffset" xOffset="0" yOffset="0" scalingFactor="1" />
  <section zOrder="1" zPosition="($CryostatLength/2)" xOffset="0" yOffset="0" scalingFactor="1" />
</xtru>
EOF
}
close(CRYOSTAT);
