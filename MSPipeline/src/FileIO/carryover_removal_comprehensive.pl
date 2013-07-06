#!/usr/bin/perl


$ns = 1;
%ordertosample = ();
%sampletoorder = ();
%list = ();
%full_list = ();
%n = ();


if($ARGV[1] eq "S") { #S: simple, Prospector output

	open(file, $ARGV[0]) or die "Can't open file: $!\n";

	while( <file> ){
		chomp;
  		@l = split /\t/, $_;
		$sample = $l[1]; #Ip ID
		$sbait{$l[1]} = $l[0]; #Bait ID
		$acc = $l[4];
		$pepcount = $l[5];
		$list{$sample}{$acc} = $pepcount;
		$withprev{$sample} = "Yes";
		++$n{$acc};
		$full_list{$sample}{$acc} = $_;
		if( !$sampletoorder{$sample} ){ 
			$ordertosample{$ns} = $sample;
			$sampletoorder{$sample} = $ns;
			$ns++;
		}
	}
}

elsif($ARGV[1] eq "F") { ##F: full, Mastertable output
	
	#open(file, "test.csv") or die "Can't open file: $!\n";
	open(file, $ARGV[0]) or die "Can't open file: $!\n";

	while( <file> ){
		chomp;
		@l = split /\t/, $_;
		$sample = $l[2]; #Ip ID
		$acc = $l[22];
		$pepcount = $l[23];
		$list{$sample}{$acc} = $pepcount;
		$sbait{$sample} = $l[7];
		++$n{$acc};
		$name{$acc} = $l[7];
		$full_list{$sample}{$acc} = $_;
		if( !$sampletoorder{$sample} ){ 
			$ordertosample{$ns} = $sample;
			$withprev{$sample} = $l[15];
			$sampletoorder{$sample} = $ns;
			$ns++;
		}
	}
}


for ( $j = 1; $j < $ns; $j++ ){
	$sample = $ordertosample{$j};
	$bait = $sbait{$sample};
	foreach $acc ( keys %{$list{$sample}} ) {
		if( $list{$sample}{$acc} > 10 ) {
			for ( $i = $j + 1; $i <= $j + 4 ; ++$i) {
				$nextsample = $ordertosample{$i};
				$nextbait = $sbait{$nextsample};
				if( $withprev{$nextsample} eq "Yes" ) {
					if ( $list{$nextsample}{$acc} > 0 && $list{$nextsample}{$acc} < $list{$sample}{$acc} / 2 && $n{$acc}<$ns/3 && $bait ne $nextbait ) {
					#if ( $list{$nextsample}{$acc} > 0 && $list{$nextsample}{$acc} < $list{$sample}{$acc} / 2 ) {
						$rem{$nextsample}{$acc} = $sample;
					}
				}
				else { 
					last; 
				};
			}
		}
	}
}

#print "Carryover\tfrom sample\tfrom bait\tto sample\tto bait\thit\n";	
for ( $i = 1; $i < $ns; $i++ ){
	$sample = $ordertosample{$i};
        foreach $acc ( keys %{$list{$sample}} ) {
		#print "Accession $acc $n{$acc}\n";
		if( !$rem{$sample}{$acc} ) {
			print "$full_list{$sample}{$acc}\n";
		}
		else{
			#print "removed $rem{$sample}{$acc}\t$sbait{$rem{$sample}{$acc}}\t$sbait{$sample}\t$full_list{$sample}{$acc}\n";
		} 
	}
}
