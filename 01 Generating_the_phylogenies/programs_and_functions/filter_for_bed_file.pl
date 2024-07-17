my @file_list = `ls *_passed.vcf`; # choose the files with "caveman_c.annot.vcf" in the title
                                               
print "@file_list \n";
                                               
foreach my $file (@file_list) {
                                               
        # open file
        print "$file";
        chomp $file;
        open (my $fh, '<', $file) or die "Can't open $file, $! \n";
        print "opened file: $file \n";
                                               
    # make an array with all the pass hits
    my @filtered = ();
    while (my $line = <$fh>) {
        if ($line =~ (/PASS/)) {
                        push @filtered, $line;
         }
    }
                                               
        # put all these hits in an output file
        my $file_out = "$file"."_for_bed_file"; # make output file with a name that links it to the input file
        print "$file_out";
        open (my $fh_out, '>', "$file_out") or die "Can't open $file_out\n";
        print $fh_out "$_" for @filtered;
        close $fh_out;
}

