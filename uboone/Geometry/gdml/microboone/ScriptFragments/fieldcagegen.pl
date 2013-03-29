#!/usr/bin/perl

$halfwidth=126;
$space=0;
$interval=1;
$i=1;

while ( $space < $halfwidth ) {
 $xcoord=$space+2;
 print '  <physvol>'."\n";
 print '   <volumeref ref="volFieldCageLoop"/>'."\n";
 print '   <position name="posFieldCageLoopPlus'."$i".'" unit="cm" x="'."$xcoord".'" y="0" z="0"/>'."\n";
 print '  </physvol>'."\n";
 print '  <physvol>'."\n";
 print '   <volumeref ref="volFieldCageLoop"/>'."\n";
 print '   <position name="posFCLoopMinus'."$i".'" unit="cm" x="-'."$xcoord".'" y="0" z="0"/>'."\n";
 print '  </physvol>'."\n";
 $space+=4*$interval;
 $i++;
}
