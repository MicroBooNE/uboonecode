#!/usr/bin/perl

$count=0;

while ( $count < 10 ) {
  $margin=$count*1.5;
  print '  <physvol> '."\n";
  print '    <volumeref ref="volElectroSetPiece"/> '."\n";
  print '    <position name="posElectroSetPieceP'."$count".'" unit="cm" x="0" y="0" z="'."$margin".'" /> '."\n";
  print '  </physvol> '."\n";
  if ( $count != 0 ) {
   print '  <physvol> '."\n";
   print '    <volumeref ref="volElectroSetPiece"/> '."\n";
   print '    <position name="posElectroSetPieceN'."$count".'" unit="cm" x="0" y="0" z="'."-$margin".'" /> '."\n";
   print '  </physvol> '."\n";
  }
  $count++;
}
