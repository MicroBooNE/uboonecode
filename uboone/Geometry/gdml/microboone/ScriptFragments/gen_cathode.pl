#!/usr/bin/perl

$i=-4;
while ( $i < 5 ) {
    $count=$i+5;
    print '  <physvol>'."\n";
    print '   <volumeref ref="volCathodeFrameVIn"/>'."\n";
    print '   <position name="posCathodeFrameVInTop'."$count".'" unit="cm" x="0" y="0.5*(kCathodeFrameWidth+kCathodeFrameVInHeight)" z="('."$i".'/5)*(0.5*kCathodeLength-kCathodeFrameWidth)"/>'."\n";
    print '  </physvol>'."\n";
    print '  <physvol>'."\n";
    print '   <volumeref ref="volCathodeFrameVIn"/>'."\n";
    print '   <position name="posCathodeFrameVInBot'."$count".'" unit="cm" x="0" y="-0.5*(kCathodeFrameWidth+kCathodeFrameVInHeight)" z="('."$i".'/5)*(0.5*kCathodeLength-kCathodeFrameWidth)"/>'."\n";
    print '  </physvol>'."\n";
    $i++;
}
