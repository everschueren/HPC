#!/usr/bin/perl                               

$pipeline ="/Users/everschueren/Projects/HPCKrogan/Scripts/Pipeline/";

#ADDED MS MACHINE column into the Master Table
# print $ARGV[0];

open(infofile, $pipeline."Gene_id_Gene_name_02_11_2011.txt") or die "Can't open file: $!\n";
while ( <infofile> ) {
  chomp;
  @l = split(/\s+/);
  $id = $l[0];
  $name = $l[1];
  $infoarray{$name} = $id; #index (gene) stores gene ID.
  $namearray{$id} = $name;
  for ($i = 2; $i < $#l+1; ++$i){
        $synname = $l[$i];
        $simarray{$ynname} = $id; #index (gene) stores gene ID.
  }  
}
close(infofile);

open(f, $pipeline."uniprot_to_entrez_human.txt") or  die "Can't open file: $!\n";

while ( <f> ) {
	chomp;
	@l = split /\t/, $_;
	$uniprot{$l[0]} = $l[2];
}
close(f);

#collapse bait names

%collapsed = ();

open(FILE,$pipeline."collapse.txt") or die "Can't open file: $!\n";
while ( <FILE> ) {
    chomp;
    @l = split(/\t/);
    $collapsed{$l[0]} = $l[1];
}
close(FILE);

#reads in Prospector input file

if($ARGV[1] eq "F"){ 
	open(f, $ARGV[0]) or  die "Can't open file: $!\n";

	while ( <f> ) {
		chomp;
		@l = split /\t/, $_;
		$l[8] =~ s/^\s+|\s+$//g; #removes white space from uniprot id
		if( $l[9] eq "Yes" ) { $l[8] = "$l[8]_$l[10]"; } #redeine bait name for co-transfected baits
		if( $collapsed{$l[8]} ) { $l[8] = $collapsed{$l[8]}; } #collapse the names of baits if specified

		$uniprotname{$l[22]} = $l[29];
		$preyorganism{$l[22]} = $l[4];
		$baitorganism{$l[8]} = $l[3];
		$proteinmw{$l[22]} = $l[27];
		#$expid{$l[8]} = $l[2]; #not one-to-one
		push @{ $exp{$l[8]}{$l[22]} }, $l[2];   
		push @{ $dataset{$l[8]}{$l[22]} }, $l[16];      
		$baitname{$l[8]} = $l[7];
	}
	close(f);
}
elsif($ARGV[1] eq "S"){
        open(f, $ARGV[0]) or  die "Can't open file: $!\n";

        while ( <f> ) {
                chomp;
                @l = split /\t/, $_;
		$s = $l[1];
		$b = $l[0];
                #print "$l[9]\t$l[10]\n";
                $b =~ s/^\s+|\s+$//g; #removes white space from uniprot id
                if( $collapsed{$b} ) { $b = $collapsed{$b}; } #collapse the names of baits if specified

                $uniprotname{$l[4]} = $l[11];
                $proteinmw{$l[4]} = $l[9];
                push @{ $exp{$b}{$l[4]} }, $s;
        	$baitname{$b} = $b;
	}
        close(f);
}

#reads Saint Output file: *_list

if ( -e "$ARGV[0]_out_saint_prob" ) { 
    # print "Reading Saint scores\n";

open(f, "$ARGV[0]_out_saint_prob");

#$fl = <f>; #first line has bait info
#	Noinfection	WT	3A	3B	3C
@fl = split /\t/, <f> ;  #first line has bait info
chomp (@fl);

while ( <f> ) {
        chomp;
        @l = split /\t/, $_;
	$prey = $l[0];
	for ($i = 1; $i < $#fl + 1; ++$i){
		$bait = $fl[$i];
		$score = $l[$i];
		$saint{$bait}{$prey} = $score;
	}
}
close(f);
}
else { 
    #print "Skipping SAINT\n";
    }

#reads CompPASS output file: *_compass

if ( -e "$ARGV[0]_compass" ) { 
    #print "Reading compPASS scores\n";

open(f, "$ARGV[0]_compass");

while ( <f> ) {
        chomp;
        @l = split /\t/, $_;
        $bait = $l[0];
        $prey = $l[1];
        $score = $l[2];
        $compass{$bait}{$prey} = $score;
}
close(f);
}
else { #print "Skipping CompPASS\n";
}



#reads MiST file with HIV training : *_mist

if ( -e "$ARGV[0]_notraining_mist.out" ) { 
#print "Reading MiST (HIV param) scores\n";

open(f, "$ARGV[0]_notraining_mist.out");
$fl = <f>; #skip first line

while ( <f> ) {
        chomp;
        @l = split /\t/, $_;
        $bait = $l[0];
        $prey = $l[1];
        $score = $l[2];
        $mist{$bait}{$prey} = $score;
}
close(f);

}
else { 
    #print "Skipping MiST (HIV param) \n";
    }
#reads MiST file with new training: *_mist

if ( -e "$ARGV[0]_training_mist.out" ) { 
#print "Reading MiST (new param) scores\n";

open(f, "$ARGV[0]_training_mist.out");
$fl = <f>; #skip first line

while ( <f> ) {
        chomp;
        @l = split /\t/, $_;
        $bait = $l[0];
        $prey = $l[1];
        $mistTR{$bait}{$prey} = $l[2];

}
close(f);
}
else { 
    #print "Skipping MiST (new param) \n";
    }

print "SampleID(s)\tDataset(s)\tBAIT\tBAITname\tBAITorganism\tPREY\tcomPASS\tMiST (HIV param)\tMiST (new param)\tUniPROT Description\tPrey Mol Weight\tPrey Organism\tEntrezID\tGeneSymbol\n";
for $b (keys %compass){  #bait
	foreach $p (keys %{$compass{$b}} ){  #prey
	$geneid = $uniprot{$p};   #print four scores and Entrez ID
        	foreach $i ( 0 .. $#{ $exp{$b}{$p} } ) {
			print "$exp{$b}{$p}[$i], ";
		}
		print "\t";
        	foreach $i ( 0 .. $#{ $dataset{$b}{$p} } ) {
			print "$dataset{$b}{$p}[$i], ";
		}
		print "\t";
		print "$b\t$baitname{$b}\t$baitorganism{$b}\t$p\t$compass{$b}{$p}\t$mist{$b}{$p}\t$mistTR{$b}{$p}\t$uniprotname{$p}\t$proteinmw{$p}\t$preyorganism{$p}\t$geneid\t$namearray{$geneid}\n";
	}
}
