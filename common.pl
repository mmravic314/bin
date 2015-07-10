sub downloadPDB {
  my $id = shift;
  my $asu = shift;
  my $ofile = shift;
  chomp($id);
  my $path;

  $id = uc($id);
  my $ml = lc(substr($id, 1, 2));

  if ($asu == 0) { $path = "data/biounit/coordinates/divided/$ml/" . lc($id) . ".pdb1.gz"; }
  else { $path = "data/structures/divided/pdb/$ml/pdb". lc($id) . ".ent.gz"; }

  if (-e "/library/databases/pdb.org/$path") {
    `gunzip -c /library/databases/pdb.org/$path > $ofile`;
  } else {
    `wget 'ftp://ftp.wwpdb.org/pub/pdb/$path' -O $ofile.gz`;
    if (-e "$ofile.gz") {
      system("cp $ofile.gz /library/databases/pdb.org/$path"); # deposit locally for later
      `gunzip -f $ofile.gz`;
    } else {
      die "error: could not find PDB file with id $id\n";
    }
  }
}

sub getCA {
  my $res = shift;

  my $ca = PDB::getAtomInRes($res, "CA", -1);
  if ($ca eq -1) {
    if ($res->{resname} =~ /(GLY|PRO|ALA|VAL|LEU|ILE|MET|CYS|PHE|TYR|TRP|HIS|HSD|HSE|HSC|HSP|LYS|ARG|GLN|ASN|GLU|ASP|SER|THR)/) {
      my ($x, $y, $z) = PDB::centerOfMass($res->{atom});
      my $ai = PDB::_newAtom(0, 0, "CA", $x, $y, $z, 0.0, 1, 0.0);
      return $ai;
    } else {
      GENERAL::warning("Failed to find CA atom for residue " . PDB::resStr($res));
      return -1;
    }
  }
  return $ca;
}

sub myBiounitPath {
  my $pdbf = shift;
  my $dirf = shift;

  return -1 if (! -d "/media/backup1/mypdb/biounits/"); # in case we are running on another machine or the device is not attached
  if (defined($dirf) && $dirf) {
    return "/media/backup1/mypdb/biounits/" . (substr(lc($pdbf), 1, 2)) . "/";
  } else {
    return "/media/backup1/mypdb/biounits/" . (substr(lc($pdbf), 1, 2)) . "/$pdbf";
  }
}

sub myDepositBiounit {
  my $spdbf = shift;
  my $dpdbf = shift;

  my $dir = myBiounitPath($dpdbf, 1);
  return if ($dir eq -1);
  GENERAL::cmkdir($dir) if (! -d $dir);
  system("cp $spdbf " . myBiounitPath($dpdbf));
}


# Generates the biological unit, given the assymetric unit, preserving only chains within the
# cutoff Ca-Ca distance from any chain in the assymetric unit.
sub generateBioUnitWithin {
  my $pdbf = shift;
  my $dcut = shift;

  my $pdb = PDB::new($pdbf);
  # Ca atoms of the asymmetric unit
  my @asuca;
  foreach my $chain (@{$pdb->{chain}}) {
    foreach my $res (@{$chain->{res}}) {
      my $ca = getCA($res); next if ($ca eq -1);
      push(@asuca, $ca);
    }
  }

  # create bio unit by apply symmetry matrices to assymmetric unit chains
  my $pdbn = PDB::new();
  my $biospec = `~/bin/mgrep.pl '(REMARK\\s+350\\s+APPLY THE FOLLOWING|REMARK\\s+350\\s+ BIOMT)' $pdbf`;
  $biospec = GENERAL::Trim($biospec);
  if ($biospec eq "") {
    return $pdb;
  }
  my @biospec = split("\n", $biospec);
  while (scalar(@biospec) > 0) {
    my $bsl = shift @biospec;
    next if ($bsl !~ /APPLY THE FOLLOWING TO CHAINS?:(.+)$/);
    my $ch = GENERAL::Trim($1);
    my @ch = split(",", $ch); @ch = GENERAL::Trim(@ch);
    if ($ch eq "") {
      GENERAL::warning("Could not parse out chain names from line '$bsl' of header for $pdbf. Will assume application to all chains!");
      foreach my $c (@{$pdb->{chain}}) { push(@ch, $c->{id}); }
    }

    my %M;
    while (scalar(@biospec) > 0) {
      last if ($biospec[0] !~ /BIOMT(\d)\s+(\d+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s+(-?[\d\.]+)\s*$/);
      my @row = ($3, $4, $5, $6);
      my $mk = $2; # matrix number
      my $rwo = $1; # row number, as provided in the matrix BIOMT? line
      my $rwe = defined($M{$mk}) ? scalar(@{$M{$mk}})+1 : 1; # expected row number
      if ($rwo != $rwe) { GENERAL::warning("Expected row $rwe, but got row $rwo... Assumming $rwe...\n"); }
      $M{$mk}->[$rwe-1] = \@row;
      shift @biospec;
    }
    # apply transformations
    foreach my $mk (keys(%M)) {
      GENERAL::assert(scalar(@{$M{$mk}}) == 3, "Matrix '$mk' parsed incorrectly from header of $pdbf");
      for (my $kk = 0; $kk < scalar(@{$M{$mk}}); $kk++) {
        GENERAL::assert(scalar(@{$M{$mk}->[$kk]}) == 4, "Matrix '$mk', row " . ($kk+1) . " has improper number of elements!");
      }
      foreach my $cid (@ch) {
        my $och = $pdb->getChain($cid, 0);
        if ($och eq 0) { print "Skipping chain $cid, not present after reading structure...\n"; next; }

        # make a separate copy of the chain
        my $pdbc = PDB::new();
        my $nch = $pdbc->_newChain($cid);
        PDB::copyChain($pdb->getChain($cid), $nch);
        my @nats = $pdbc->conAtoms();

        # apply transform
        foreach my $na (@nats) {
          my $x = $na->{xcoor}; my $y = $na->{ycoor}; my $z = $na->{zcoor};
          $na->{xcoor} = $x*($M{$mk}->[0]->[0]) + $y*($M{$mk}->[0]->[1]) + $z*($M{$mk}->[0]->[2]) + ($M{$mk}->[0]->[3]);
          $na->{ycoor} = $x*($M{$mk}->[1]->[0]) + $y*($M{$mk}->[1]->[1]) + $z*($M{$mk}->[1]->[2]) + ($M{$mk}->[1]->[3]);
          $na->{zcoor} = $x*($M{$mk}->[2]->[0]) + $y*($M{$mk}->[2]->[1]) + $z*($M{$mk}->[2]->[2]) + ($M{$mk}->[2]->[3]);
        }

        # check that the chain is within cutoff distance of ASU after transformation -- build a box containing all CA atoms of the new chain
        my $xmi = 10**10; my $xma = -10**10; my $ymi = $xmi; my $yma = $xma; my $zmi = $xmi; my $zma = $xma;
        my $nca = 0;
        foreach my $res (@{$nch->{res}}) {
          my $ca = getCA($res); next if ($ca eq -1);
          $xmi = GENERAL::min($xmi, $ca->{xcoor}); $xma = GENERAL::max($xma, $ca->{xcoor});
          $ymi = GENERAL::min($ymi, $ca->{ycoor}); $yma = GENERAL::max($yma, $ca->{ycoor});
          $zmi = GENERAL::min($zmi, $ca->{zcoor}); $zma = GENERAL::max($zma, $ca->{zcoor});
          $nca++;
        }
        next if ($nca == 0); # skip chains with no CA atoms
        my $cutf = 0; # see if distance from any CA atom of the ASU to the box around the transformed chain is below cutoff
        foreach my $ca (@asuca) {
          if (($ca->{xcoor} > $xmi-$dcut) && ($ca->{xcoor} < $xma+$dcut) && ($ca->{ycoor} > $ymi-$dcut) && ($ca->{ycoor} < $yma+$dcut) && ($ca->{zcoor} > $zmi-$dcut) && ($ca->{zcoor} < $zma+$dcut)) {
            $cutf = 1;
            last;
          }
        }
        if (($cutf == 0) && ($M{$mk}->[0]->[0] == 1) && ($M{$mk}->[1]->[1] == 1) && ($M{$mk}->[2]->[2] == 1) && ($M{$mk}->[0]->[3] == 0) && ($M{$mk}->[1]->[3] == 0) && ($M{$mk}->[2]->[3] == 0)) {
          die("weird, a chain in the ASU did not make it after cutting ($pdbf)...\n");
        }
        next if ($cutf == 0); # skip chain if it is out of range

        # if within cutoff, copy chain to bio unit
        if ($cutf) {
          my $nbuc = $pdbn->_newChain($cid);
          PDB::copyChain($nch, $nbuc);
        }
      }
    }
  }
  return $pdbn;
}


# finds the smallest sphere that contains all atoms in the given array
sub getCircumscribingSphere {
  my $atoms = shift;

  if (scalar(@$atoms) == 0) {
    return (0, 0, 0, 0);
  }

  my $cx = 0; my $cy = 0; my $cz = 0;
  foreach my $a (@$atoms) {
    $cx += $a->{xcoor};
    $cy += $a->{ycoor};
    $cz += $a->{zcoor};
  }
  $cx /= scalar(@$atoms);
  $cy /= scalar(@$atoms);
  $cz /= scalar(@$atoms);

  my $R = 0;
  foreach my $a (@$atoms) {
    $R = GENERAL::max($R, sqrt(($a->{xcoor} - $cx)**2 + ($a->{ycoor} - $cy)**2 + ($a->{zcoor} - $cz)**2));
  }
  return ($cx, $cy, $cz, $R);
}

1;
