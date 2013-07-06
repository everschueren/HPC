#!/usr/bin/perl


if ($ARGV[2] eq "F") {

#open(FILE, "new_additions_hip44-56_key.txt") or  die "Can't open file: $!\n";
	open(FILE, $ARGV[0]) or  die "Can't open file: $!\n";

	$fl = <FILE>;

	while ( <FILE> ) {
		chomp;
		@l = split /\t/, $_;
		$id = "$l[2]_$l[7]"; #column 19 to be filled
		$store{$l[2]} = "$_$id";
	}
	close(FILE);

#open(FILE, "new_additions_hip44-56_data.txt") or  die "Can't open file: $!\n";
	open(FILE, $ARGV[1]) or  die "Can't open file: $!\n";

	print "\t#\n";
	while ( <FILE> ) {
		chomp;
		@l = split /\t/, $_;

		print "$store{$l[0]}\t$l[1]\t$l[2]\t$l[3]\t$l[4]\t$l[5]\t$l[6]\t$l[7]\t$l[8]\t$l[9]\t$l[10]\n";
	}
	close(FILE);
}

elsif ($ARGV[2] eq "S"){

#open(FILE, "key.txt") or  die "Can't open file: $!\n";
        open(FILE, $ARGV[0]) or  die "Can't open file: $!\n";

        $fl = <FILE>;

        while ( <FILE> ) {
                chomp;
                @l = split /\t/, $_;
                $store{$l[0]} = $l[1];
        }
        close(FILE);

#open(FILE, "new_additions_hip44-56_data.txt") or  die "Can't open file: $!\n";
        open(FILE, $ARGV[1]) or  die "Can't open file: $!\n";
	print "\t#\n";

        while ( <FILE> ) {
                chomp;
                @l = split /\t/, $_;

                print "$store{$l[0]}\t$_\n";
        }
        close(FILE);

}

