#!/usr/bin/perl

print "usage: SimplifyEvidence <evidence.txt file> <files.txt>\n";
print "files.txt is a tab-delimited text file providing sample names for each raw file name\n";

open (infile, $ARGV[1]);
while ($line = <infile>) {
    chomp ($line);
    @pieces = split /\t/, $line;
    $rawfile = $pieces[0];
    $sample = $pieces[1];
    $samplehash{$rawfile} = $sample;
    print "$rawfile\t$samplehash{$rawfile}\n";
}
close (infile);

open (infile, $ARGV[2]);
$numfilters = 0;
while ($line = <infile>) {
    chomp ($line);
    @pieces = split /\t/, $line;
    $filtercolumn[$numfilters] = $pieces[0];
    $filtercondition[$numfilters] = $pieces[1];
    print "$filtercolumn[$numfilters]\t$filtercondition[$numfilters]\n";
    $numfilters++;
}

open (infile, $ARGV[0]);
$line = <infile>;
@pieces = split /\t/, $line;
for ($i=0; $i<=$#pieces; $i++) {
    $columns{$pieces[$i]} = $i;
}

if ($columns{"Raw file"}) {
    $rawfilecolumn = $columns{"Raw file"};
} elsif ($columns{"Raw File"}) {
    $rawfilecolumn = $columns{"Raw File"};
}
if ($columns{"Mod. peptide ID"}) {
    $modpepidcolumn = $columns{"Mod. peptide ID"};
} elsif ($columns{"Mod. Peptide ID"}) {
    $modpepidcolumn = $columns{"Mod. Peptide ID"};
}
if ($columns{"Protein Descriptions"}) {
    $descriptioncolumn = $columns {"Protein Descriptions"};
} elsif ($columns{"Fasta headers"}) {
    $descriptioncolumn = $columns{"Fasta headers"};
}
$ratiocolumn = $columns{"Ratio H/L"};
if ($columns{"Modified sequence"}) {
    $modseqcolumn = $columns{"Modified sequence"};	   
} elsif ($columns{"Modified Sequence"}) {
    $modseqcolumn = $columns{"Modified Sequence"};
}
$proteinscolumn = $columns{"Proteins"};

while ($line = <infile>) {
    @pieces = split /\t/, $line;
    $modpepid = $pieces[$modpepidcolumn];
    $sample = $samplehash{$pieces[$rawfilecolumn]};
    $description{$modpepid} = $pieces[$descriptioncolumn];
    $pepid = $modpepid . "." . $sample;
    $proteins{$modpepid} = $pieces[$proteinscolumn];
    if ($pieces[$modseqcolumn] =~ /\S+/) {
	$modseq{$modpepid} = $pieces[$modseqcolumn];
    }

    print "$modpepid\t$sample\t$description{$modpepid}\n";

    if ($pieces[$ratiocolumn] =~ /\d+/) {
	if ($numratios{$pepid}) {
	    $logratios{$pepid}[$numratios{$pepid}] = log($pieces[$ratiocolumn]);
	    $numratios{$pepid}++;
	} else {
	    $logratios{$pepid}[0] = log($pieces[$ratiocolumn]);
	    $numratios{$pepid} = 1;
	}
    } elsif ($pieces[$ratiocolumn] =~ /\S+/) {
	if ($numnanvalues{$pepid}) {
	    $numnanvalues{$pepid}++;
	} else {
	    $numnanvalues{$pepid} = 1;
	}
	print "NaN $numnanvalues{$pepid}\n";

	if (!($numratios{$pepid})) {
	    $numratios{$pepid} = 0;
	}
    }
}
close (infile);

$outfilename = ">" . substr $ARGV[0], 0, length($ARGV[0])-4;
$outfilename = $outfilename . "-simple.txt";
open (outfile, $outfilename);

print outfile "Modified Peptide Sequence\tProteins\t$Protein Description\tSample\tAvg Log Ratio\tNorm Avg Log Ratio\tZ-score\tNum Ratios\tNum NaN Values\n";

foreach $key (keys %numratios) {
    $key =~ /(\S+)\.(\S+)/;
    $modpepid = $1;
    $sample = $2;
    
    if ($numratios{$key} > 0) {
	$sum = 0;
	for ($i=0; $i<$numratios{$key}; $i++) {
	    $sum += $logratios{$key}[$i];
	}
	$avg = $sum / $numratios{$key};
	
	if ($samplenumratios{$sample}) {
	    $sampleratios{$sample}[$samplenumratios{$sample}] = $avg;
	    $samplenumratios{$sample}++;
	} else {
	    $sampleratios{$sample}[0] = $avg;
	    $samplenumratios{$sample} = 1;
	}
    }
}

foreach $key (keys %samplenumratios) {
    $sum = 0;
    for ($i=0; $i<$samplenumratios{$key}; $i++) {
	$sum += $sampleratios{$key}[$i];
    }
    $sampleavg{$key} = $sum / $samplenumratios{$key};
    
    $sumsquares = 0;
    for ($i=0; $i<$samplenumratios{$key}; $i++) {
	$sumsquares += ($sampleratios{$key}[$i] - $sampleavg{$key}) * ($sampleratios{$key}[$i] - $sampleavg{$key});
    }
    $samplestdev{$key} = sqrt($sumsquares / ($samplenumratios{$key} - 1));
    
    print "$key $sampleavg{$key} $samplestdev{$key}\n";
}


foreach $key (keys %numratios) {
    $key =~ /(\S+)\.(\S+)/;
    $modpepid = $1;
    $sample = $2;
    
    print outfile "$modseq{$modpepid}\t$proteins{$modpepid}\t$description{$modpepid}\t$sample\t";
    
    if ($numratios{$key} > 0) {
	$sum = 0;
	for ($i=0; $i<$numratios{$key}; $i++) {
	    $sum += $logratios{$key}[$i];
	}
	$avg = $sum / $numratios{$key};
	$normavg = $avg - $sampleavg{$sample};
	$zscore = $normavg / $samplestdev{$sample};
	print outfile "$avg\t$normavg\t$zscore\t$numratios{$key}\t";
    } else {
	print outfile "NA\t0\t0\t0\t";
    }

    if ($numnanvalues{$key}) {
	print outfile "$numnanvalues{$key}\n";
    } else {
	print outfile "0\n";
    }
    


}
close (infile);
close (outfile);
	    
	    


	
	
