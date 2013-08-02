#!/usr/bin/perl

#Version 2.0--Sept 12, 2012
#Disregards certain samples specified in "samples_to_remove.txt" excluded due to low expression, poor quality etc.
#Excludes specificity between certain baits specified in "specificity_exclusions.txt" file, GP160-GP120, mutant version of the bait, experimental conditions such as MG132 etc.
#contact: natali.gulbahce@gmail.com

%rem_smp = ();
$data = $ARGV[0];
$format = $ARGV[1];
$remove = $ARGV[2];
$collapse = $ARGV[3];
$exclusions = $ARGV[4];

# print $remove;
open(FILE,$remove) or die "Can't open file: $!\n";
while ( <FILE> ) {
    chomp;
    @l = split(/\t/);
    $rem_smp{$l[0]} = 1;
}
close(FILE);

#collapse bait names

%collapsed = ();

open(FILE,$collapse) or die "Can't open file: $!\n";
while ( <FILE> ) {
    chomp;
    @l = split(/\t/);
    $collapsed{$l[0]} = $l[1];
}
close(FILE);

%seen = ();
%seen_bait = ();
%numuniq = ();

#reads in the IP-MS file

if($format eq "F"){
	open(msfile, $data) or die "Can't open file: $ARGV[0]\n";
	$fl = <msfile>;

	while ( <msfile> ) {
		chomp;
		@l = split(/\t/);
		

		if(!$rem_smp{$l[2]} ){ #if the sample is of good quality continue
			$l[8] =~ s/^\s+|\s+$//g; #removes white space from uniprot id
			#$l[2] =~ s/^\s+|\s+$//g; #removes white space from sample id
			#$l[22] =~ s/^\s+|\s+$//g; #removes white space from prey uniprot
			if( $l[9] eq "Yes" ) { $l[8] = "$l[8]_$l[10]"; } #redefine bait name for co-transfected baits
			if( $collapsed{$l[8]} ) { $l[8] = $collapsed{$l[8]}; } #collapse the names of baits if specified
	 
	 		if( $l[21] > 0 || $l[21] eq "" ) { #require uniq peptides
				$length{$l[22]} = int($l[27]/110) ; #stores protein length in prey uniprot id
				if(!$seen{$l[2]}) {
					push @samplelist, $l[2];
					$seen{$l[2]} = 1;
					$seen_bait{$l[8]} = 1;	#if bait is not seen than specificity exclusion for it will fail below
				}
				$numuniq{$l[22]}{$l[2]} = $l[23];  #stores peptide counts per uniprot and sample	
				$idbait{$l[2]} = $l[8];  # bait id in uniprot
				$namebait{$l[2]} = $l[7]; # bait id in uniprot name
			}
		}
	}
}
elsif($format eq "S"){
	open(msfile, $data) or die "Can't open file: $!\n";
	$fl = <msfile>;

	while ( <msfile> ) {
		chomp;
		@l = split(/\t/);
		$b = $l[0];
		$s = $l[1];


		if(!$rem_smp{$s} ){ #if the sample is of good quality continue

			$b =~ s/^\s+|\s+$//g; #removes white space from uniprot id
			
			if( $collapsed{$b} ) { $b = $collapsed{$b}; } #collapse the names of baits if specified
	 
			if( $l[3] > 0 || $l[3] eq "" ) { #require uniq peptides
				$length{$l[4]} = int($l[9]/110) ; #stores protein length in prey uniprot id
				if(!$seen{$s}) {
					push @samplelist, $s;
					$seen{$s} = 1;
					$seen_bait{$b} = 1;	#if bait is not seen than specificity exclusion for it will fail below
				}
				$numuniq{$l[4]}{$s} = $l[5];  #stores peptide counts per uniprot and sample	
				$idbait{$s} = $b;  # bait id in uniprot
				$namebait{$b} = $l[11]; # bait id in uniprot name
			}
		}else{
			#print "REMOVING\t".$s
		}
	}
}



open(FILE, $exclusions) or die "Can't open file: $!\n";
while(<FILE>){
        chomp;
        @l = split(/\t/);
        $spec = "";
        if ($seen_bait{$l[0]}) {
                @e = split('\|',$l[1]);
                for ($i = 0; $i <= $#e; $i++){ 
                        if ($seen_bait{$e[$i]} ) { $spec .= "$e[$i]"; $flag = 1;} 
                        if ($i == $#e ) { $flag = 0}
                        if ($flag) { $spec .= "|"; $flag = 0;}
                }
                if($spec =~ /\|$/) { chop($spec); } 
                $exc_bait{$l[0]} = $spec;
                #print "$l[0] $spec\n";
        }
}

#first line
print "a\tb\tc\tIP";
for ($i = 0; $i <= $#samplelist; $i++) {
	print "\t$samplelist[$i]";
}
print "\n";

#second line
print "a\tb\tc\tBait";
for ($i = 0; $i <= $#samplelist; $i++) {
	$id = $samplelist[$i];
        print "\t$idbait{$id}";
}
print "\n";

#third line
print "Preys\tPepAtlas\tLength\tPreyType\\BaitCov";
for ($i = 0; $i <= $#samplelist; $i++) {
        $id = $samplelist[$i];
	$bait = $idbait{$id};
	if($exc_bait{$bait}) { print "\t$exc_bait{$bait}"; }
	else { print "\t$bait"; }
}
print "\n";

for $prot (keys %numuniq){
	print "$prot\t0\t$length{$prot}\tN";
	for ($i = 0; $i <= $#samplelist; $i++) {
		$id = $samplelist[$i];
		if($numuniq{$prot}{$id}) { print "\t$numuniq{$prot}{$id}";}
		else { print "\t0"; }
	}
	print "\n";
}



