#!/usr/bin/perl

@EXPORT = qw"PDBline";

sub PDBline {
    my $line = $_[0];
    my $field = $_[1];
    my $set = 1 if defined $_[2];
    my $value = $_[2] if $set;

    return $line unless( ($line =~ /^HETATM/) or ($line =~ /^ATOM/));#if not an atom, return line

    $START = 0;
    $END   = 1;
    %columns = (
        atom        => [1,6],
        atomserial  => [7,11],
        spacer1     => [12,12],
        atomname    => [13,16],
	altLoc      => [17,17],
        resname     => [18,20],
        spacer2     => [21,21],
        chainid     => [22,22],
        resseq      => [23,26],
        icode       => [27,27],
        spacer3     => [28,30],
        x           => [31,38],
        y           => [39,46],
        z           => [47,54],
        occupancy   => [55,60],
        tempfactor  => [61,67],#changed fromm 72 because spacer4 was uncommented
                spacer4     => [67,72],
        segid       => [73,76],
        element     => [77,78],
        totalcharge => [79,80]
        );



    return substr($line,$columns{$field}[$START]-1,$columns{$field}[$END]-$columns{$field}[$START]+1) if (!$set);

    $length = $columns{$field}[$END] - $columns{$field}[$START]+1;
    $formated_value = substr($value,0,$length);

    if (length($value) < $length) {
    $i = $length - length($value);
    while ($i != 0) {
        $formated_value = " ".$formated_value;#changed by S.Hershman 6/28/06 from $formated_value." "
        $i --;
    }
    }
    $newline = substr($line,0,$columns{$field}[$START]-1).$formated_value.substr($line,$columns{$field}[$END]);

    return $newline;


}

#         1         2         3         4         5         6         7         8
#12345678901234567890123456789012345678901234567890123456789012345678901234567890
#ATOM     12  HD22ASN     1       3.360   5.504   2.994  0.00  0.00           H
#ATOM      1  CA  ARG     1       3.970   2.846   0.000  0.00  0.00           C


%aminos = (
        "ALA" => {
        one  => "A",
        full => "Alanine",
       },
       "ARG" => {
           one  => "R",
           full => "Arginine"
       },
       "ASN" => {
           one  => "N",
           full => "Asparagine",
       },
           "ASP" => {
           one  => "D",
           full => "Aspartic acid",
       },
           "CYS" => {
           one  => "C",
           full => "Cysteine",
       },
           "GLN" => {
           one  => "Q",
           full => "Glutamine",
       },
           "GLU" => {
           one  => "E",
           full => "Glutamic acid",
       },
           "GLY" => {
           one  => "G",
           full => "Glycine",
       },
           "HIS" => {
           one  => "H",
           full => "Histidine",
       },
           "ILE" => {
           one  => "I",
           full => "Isoleucine",
       },
           "LEU" => {
           one  => "L",
           full => "Leucine",
       },
           "LYS" => {
           one  => "K",
           full => "Lysine",
       },
           "MET" => {
           one  => "M",
           full => "Methionine",
       },
           "PHE" => {
           one  => "F",
           full => "Phenylalanine",
       },
           "PRO" => {
           one  => "P",
           full => "Proline",
       },
           "SER" => {
           one  => "S",
           full => "Serine",
       },
           "THR" => {
           one  => "T",
           full => "Threonine",
       },
           "TRP" => {
           one  => "W",
           full => "Tryptophan",
       },
           "TYR" => {
           one  => "Y",
           full => "Tyrosine",
       },
           "VAL" => {
           one  => "V",
           full => "Valine",
       },
           "XAA" => {
           one  => "X",
           full => "Any or unknown amino acid",
       },
       );

