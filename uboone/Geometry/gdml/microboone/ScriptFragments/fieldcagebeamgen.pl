#!/usr/bin/perl

$i=-5;
$count=7;

while ( $i < 6 ) {
     $fraction=$i/5;
     print '  <physvol>'."\n";
     print '   <volumeref ref="volFieldCageBarFlat"/>'."\n";
     print '   <position name="posFieldCageBarTop'."$count".'" unit="cm" x="0" y="kFieldCageBeamYPos" z="'."$fraction".'*kFieldCageBeamZInt"/>'."\n";
     print '   <rotation name="rFieldCageBarTop'."$count".'" unit="cm" x="90" y="0" z="0"/>'."\n";
     print '  </physvol>'."\n";
     print '  <physvol>'."\n";
     print '   <volumeref ref="volFieldCageBarFlat"/>'."\n";
     print '   <position name="posFieldCageBarBot'."$count".'" unit="cm" x="0" y="-kFieldCageBeamYPos" z="'."$fraction".'*kFieldCageBeamZInt"/>'."\n";
     print '   <rotation name="rFieldCageBarBot'."$count".'" unit="cm" x="-90" y="0" z="0"/>'."\n";
     print '  </physvol>'."\n";
     $i++;
     $count++;
}
