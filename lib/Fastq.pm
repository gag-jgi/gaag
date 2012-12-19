package JGI::Bio::Utils::Fastq;

require Exporter;
@ISA = qw(Exporter);
@EXPORT =
  qw(determineFileType cntReads cntEstimatedReads fqReadLen splitFq splitFqRandomly medianReadLen convertToIlluminaQuals selectRandSeqs isPairedEnd separatePairedEndFq
  cntSeqs splitFaRandomly selectSeqsById filterByLen calcContigStatsMinLen separatePairedEndFa splitSingleLineFastaRandomly fastqToFastaQual createBwaIndex
  guessQualityScoreEncoding renameSeqs splitFq splitFa decodeQuals decodeSangerQuals selectRandSeqsPerl);

use strict;
use warnings;

use Bio::SeqIO;
use File::Basename;
use JGI::Bio::Utils;
use JGI::Bio::Utils::Math;
use Data::Alias;
use File::Temp qw/ tempfile tempdir /;

our $VERSION = '1.3.1';

=head1 NAME

JGI::Bio::Utils::Fastq
=head1 SYNOPSIS

use JGI::Bio::Utils::Fastq;
my $read_length = fqReadLen(in.fq);

=head2 DESCRIPTION

Utils for common operations relating to FASTA and FASTQ files

=head1 SEE ALSO

JGI::Bio::Rnnotator::ReadTrimmer
JGI::Bio::Rnnotator::ReadPreprocessor

=head1 AUTHOR

Jeff Martin, E<lt>jamartin@lbl.govE<gt>
Zhong Wang, E<lt>zhongwang@lbl.govE<gt>
Xiandong Meng, E<lt>xiandongmeng@lbl.govE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011-2012 by JGI

=cut

# Checks the given file to determine if the reads are in FASTA or FASTQ format.
# 	@param $file - the file path
# 	@returns - 'fastq', 'fasta', or 'unknown'
sub determineFileType($) {
    my ($file) = @_;

    open(IN, $file) || die("Unable to open input file: $file\n");
    my $first_line = <IN>;
    my $first_char = substr($first_line, 0, 1);
    if ($first_char eq '@') {
        return 'fastq';
    }
    elsif ($first_char eq '>') {
        return 'fasta';
    }
    else {
        return 'unknown';
    }
    close(IN);
}

##############################################
# 			         FASTQ
##############################################

# Creates a BWA index using the appropriate algorithm
sub createBwaIndex {
    my ($ref) = @_;
    defined($ref) || die("No ref defined.");
    (-e $ref) || die("Ref does not exist: $ref");
    my ($tmp_file, $tmp_dir) = fileparse($ref);
    (-w $tmp_dir) || die("Ref. directory is not writable: $tmp_dir");
    if (-s $ref > 2147483648) {    # 2 GB
        system("bwa index $ref -a bwtsw > /dev/null 2>&1") unless (-e "$ref.bwt");
    }
    else {
        system("bwa index $ref > /dev/null 2>&1") unless (-e "$ref.bwt");
    }
}

# Note: I stoll this from Joel (/home/jel/bin/fastq_guess_format)
# not checking between solexa and illumiina
# returns: uknown, sanger, or illumina
sub guessQualityScoreEncoding {
    my $file = shift;

    if (!defined $file || ($file ne '-' && !-e $file)) {
        print STDERR "$0 file.fq[.gz, .bz2]\n";
        exit;
    }

    my $fh;
    if ($file eq '-') {
        $fh = \*STDIN;
    }
    else {
        fqopen($fh, $file) or die $!;
    }

    my $type = "unknown";
    while (<$fh>) {
        <$fh>;
        <$fh>;
        $_ = <$fh>;    # fast forward to quality
        if (/[!-:]/) {
            return 'sanger';
        }
        elsif (/[K-h]/) {
            return 'illumina';
        }
    }
    return $type;
}

# Note: I stoll this from Joel (/home/jel/bin/fastq_guess_format)
sub fqopen {
    alias my ($fh, $file) = @_;

    my $result;
    my $open =
        $file =~ /[.]gz$/  ? "gzip  -cd $file |"
      : $file =~ /[.]bz2$/ ? "bzip2 -cd $file |"
      :                      $file;

    return (open($fh, $open));
}

# Converts a string of encoded quality scores to their numeric format.
# 	@param $quality_str - the string containing the quality scores
# 	@returns - an array ref to the numeric quality scores
sub decodeQuals {
    my ($quality_str, $offset) = @_;
    my @result = unpack("C*", $quality_str);
    for (my $i = 0 ; $i < @result ; $i++) {
        $result[$i] -= $offset;
    }
    return \@result;
}

# Converts a string of encoded quality scores to their numeric format.
# 	@param $quality_str - the string containing the quality scores
# 	@returns - an array ref to the numeric quality scores
sub decodeSangerQuals {
    my ($quality_str) = @_;
    decodeQuals($quality_str, 33);
}

# Converts a string of encoded quality scores to their numeric format.
# 	@param $quality_str - the string containing the quality scores
# 	@returns - an array ref to the numeric quality scores
sub convertToIlluminaQuals {
    my ($quality_str) = @_;
    decodeQuals($quality_str, 64);
}

# Splits a FASTQ file
sub splitFq {
    my ($filename, $num_out, $outdir) = @_;
    my $pe = isPairedEnd($filename);

    # set outdir if not already set
    if (!defined($outdir)) {
        my ($tmp_file, $tmp_dir) = fileparse($filename);
        $outdir = $tmp_dir;
    }
    my $base = basename($filename);

    if ($num_out == 0) {
        die("Cannot split a file into 0 jobs.");
    }

    # determine how many reads we want in each file
    my $num_reads          = cntEstimatedReads($filename);
    my $num_reads_per_file = int($num_reads / $num_out);

    # split
    my $curr_file_cnt = 1;
    my $curr_read_cnt = 0;
    my $curr_fh;
    open($curr_fh, ">$outdir/$base.$curr_file_cnt")
      || die("Can't write: $outdir/$base.$curr_file_cnt");
    open(IN, $filename) || die($!);
    my @split_files = ("$outdir/$base.$curr_file_cnt");
    while (defined(my $h1 = <IN>)
        && defined(my $s  = <IN>)
        && defined(my $h2 = <IN>)
        && defined(my $q  = <IN>))
    {
        $curr_read_cnt++;
        print $curr_fh $h1, $s, $h2, $q;
        if (   $pe
            && defined(my $h1b = <IN>)
            && defined(my $sb  = <IN>)
            && defined(my $h2b = <IN>)
            && defined(my $qb  = <IN>))
        {
            print $curr_fh $h1b, $sb, $h2b, $qb;
            $curr_read_cnt++;
        }

        if ($curr_read_cnt >= $num_reads_per_file) {
            close($curr_fh);
            $curr_file_cnt++;
            $curr_read_cnt = 0;
            open($curr_fh, ">$outdir/$base.$curr_file_cnt");
            push(@split_files, "$outdir/$base.$curr_file_cnt");
        }
    }
    close(IN);
    close($curr_fh) if ($curr_fh);

    return \@split_files;
}

# Randomly splits a FASTQ file
sub splitFqRandomly {
    my ($filename, $num_out, $pe) = @_;

    if ($num_out == 0) {
        die("Cannot split a file into 0 jobs.");
    }

    # open output files
    my @split_files = ();
    my @out_fhs     = ();
    for (my $file_cnt = 1 ; $file_cnt <= $num_out ; $file_cnt++) {
        my $split_file = "$filename.$file_cnt";
        push(@split_files, $split_file);
        my $fh;
        open($fh, ">$split_file") || die("Unable to open $split_file for writing\n");
        push(@out_fhs, $fh);
    }

    # write to several output files
    open(IN, $filename) || die("Unable to open input file: $filename\n");
    while (defined(my $h1 = <IN>)
        && defined(my $s  = <IN>)
        && defined(my $h2 = <IN>)
        && defined(my $q  = <IN>))
    {

        # randomly choose an output file
        my $out_fh = $out_fhs[ int(rand($num_out)) ];

        # write to the appropriate file (ensures every sequence is output somewhere)
        print $out_fh $h1, $s, $h2, $q;
        if (   $pe
            && defined(my $h1b = <IN>)
            && defined(my $sb  = <IN>)
            && defined(my $h2b = <IN>)
            && defined(my $qb  = <IN>))
        {
            print $out_fh $h1b, $sb, $h2b, $qb;
        }
    }
    close(IN);

    # close output files
    for my $out_fh (@out_fhs) {
        close($out_fh);
    }

    return \@split_files;
}

# Randomly select $num of reads from the given fastq file.
# 	@param $fq - the fastq file
# 	@param $num - the nubmer of reads to randomly select (paired-end counted as two)
# 	@param $outfile - the output file containing the randomly selected reads
# 	@returns - (nothing)
sub selectRandSeqs {
    my ($fq, $num, $outfile) = @_;

    # use Rob's ultra-fast sampler
    my $pe = isPairedEnd($fq);

    # divide by two since RandomlySample selects twice too many reads for PE runs
    if ($pe) {
        $num /= 2;
        $num = int($num);
    }
    my $cmd = "RandomlySample --num $num $fq > $outfile";
    system($cmd);
    unless (-e $outfile && -s $outfile > 0) {
        die("Output file from random selection is empty... using command: $cmd.");
    }
}

# Randomly select $num of reads from the given fastq file using perl.
# 	@param $fq - the fastq file
# 	@param $num - the nubmer to randomly select
# 	@param $outfile - the output file containing the randomly selected reads
# 	@returns - (nothing)
sub selectRandSeqsPerl {
    my ($fq, $num, $outfile) = @_;
    my $pe    = isPairedEnd($fq);
    my $reads = cntEstimatedReads($fq);

    # update to number of fragments
    if ($pe) {
        $reads /= 2;
        $num   /= 2;
    }
    if ($num > $reads) {
        die("Not enough reads in FASTQ to randomly select. $reads reads in FASTQ. Selecting $num\n"
        );
    }

    # generate random numbers
    my %sel;
    for (my $i = 0 ; $i < $num ; $i++) {
        while (my $r = int(rand($reads) + 1)) {
            next if ($sel{$r});
            $sel{$r} = 1;
            last;
        }
    }

    # stream through the FASTQ file
    open(IN, $fq) or die $!;
    open(OUT, ">$outfile");
    my $n = 0;
    if ($pe) {
        while ( defined(my $h1 = <IN>)
            and defined(my $s1  = <IN>)
            and defined(my $qh1 = <IN>)
            and defined(my $q1  = <IN>)
            and defined(my $h2  = <IN>)
            and defined(my $s2  = <IN>)
            and defined(my $qh2 = <IN>)
            and defined(my $q2  = <IN>))
        {
            $n++;
            print OUT $h1, $s1, $qh1, $q1, $h2, $s2, $qh2, $q2 if ($sel{$n});
        }
    }
    else {
        while ( defined(my $h1 = <IN>)
            and defined(my $s1  = <IN>)
            and defined(my $qh1 = <IN>)
            and defined(my $q1  = <IN>))
        {
            $n++;
            print OUT $h1, $s1, $qh1, $q1 if ($sel{$n});
        }
    }
    close(IN);
    close(OUT);
}

# Counts the nubmer of reads in the given fastq file
# 	@param $fq - the fastq file
# 	@returns - the number of reads
sub cntReads {
    my ($fq) = @_;
    return (num_lines($fq) / 4);
}

# Estimate the approximate number of reads in an Illumina fastq
# 	@param $fq - the fastq file
# 	@returns - the aproximate number of reads
sub cntEstimatedReads {
    my ($fq) = @_;

    die("Fastq file does not exist: $fq") unless (-e $fq);
    if (-s $fq == 0) {
        return 0;
    }
    my $result = `EstimatedReadNums $fq`;
    my @result = split(/\s+/, $result);
    my $answer = $result[4];
    defined($answer)
      || die("Estimating number of reads failed when executing: EstimatedReadNums $fq\n");
    return $answer;
}

# Randomly samples 100 reads and determines the read length based upon the random sampling. Returns -1 if
# the reads are not the same length.
# 	@param $fq - the fq file
# 	@returns - the read length
sub fqReadLen {
    my ($fq) = @_;

    # randomly select seqs
    my $tmp_dir = tempdir(CLEANUP => 1);
    my ($fh, $rand) = tempfile(DIR => $tmp_dir);
    selectRandSeqs($fq, 100, $rand);

    # check read lengths
    open(IN, $rand) || die("Unable to open fq: $rand\n");
    my $read_len;
    while ( defined(my $h1 = <IN>)
        and defined(my $s  = <IN>)
        and defined(my $h2 = <IN>)
        and defined(my $q  = <IN>))
    {
        chomp($s);
        my $curr_read_len = length($s);
        if (!defined($read_len)) {
            $read_len = $curr_read_len;
        }
        elsif ($curr_read_len != $read_len) {
            close(IN);
            return -1;
        }
    }
    close(IN);

    unlink($rand);
    return $read_len;
}

# Randomly samples 5000 reads to determine the median read length
# 	@param $fq - the fq file
# 	@returns - ($read1_median_len, $read2_median_len), second param is -1 in the single-end case
sub medianReadLen {
    my ($fq) = @_;

    my $tmp_dir = tempdir(CLEANUP => 1);
    my ($fh, $rand) = tempfile(DIR => $tmp_dir);

    my $pe = isPairedEnd($fq);
    selectRandSeqs($fq, 5_000, $rand);
    open(IN, $rand) || die("Unable to read fq: $rand\n");
    my @read1_lengths = ();
    my @read2_lengths = ();
    while (defined(my $h = <IN>)
        && defined(my $s  = <IN>)
        && defined(my $h2 = <IN>)
        && defined(my $q  = <IN>))
    {
        chomp($s);
        push(@read1_lengths, length($s));
        if (   $pe
            && defined(my $hb  = <IN>)
            && defined(my $sb  = <IN>)
            && defined(my $h2b = <IN>)
            && defined(my $qb  = <IN>))
        {
            chomp($sb);
            push(@read2_lengths, length($sb));
        }
    }
    unlink($rand);
    close(IN);

    my $read1_median_len = median(\@read1_lengths);
    my $read2_median_len = @read2_lengths > 0 ? median(\@read2_lengths) : -1;
    return ($read1_median_len, $read2_median_len);
}

# Looks at the sequence identifiers to determine if the given fastq file is single or paired end
# 	@param $fq - the fq file
# 	@returns - true if the fq is paired end, false otherwise
sub isPairedEnd {
    my ($fq) = @_;
    open(IN, $fq) or die("Unable to open fq: $fq\n");
    my $is_paired_end = 0;
    if (    defined(my $h1 = <IN>)
        and defined(my $s    = <IN>)
        and defined(my $h2   = <IN>)
        and defined(my $q    = <IN>)
        and defined(my $h1_b = <IN>)
        and defined(my $s_b  = <IN>)
        and defined(my $h2_b = <IN>)
        and defined(my $q_b  = <IN>))
    {
        chomp($h1);
        chomp($h1_b);
        my @h1_words = split(' ', $h1);
        $h1 = $h1_words[0];
        my @h1_b_words = split(' ', $h1_b);
        $h1_b = $h1_b_words[0];
        $is_paired_end = (($h1 =~ m/.+\/1$/) && (($h1_b =~ m/.+\/2$/))) ? 1 : 0;

	if (!$is_paired_end){
		my $h11 = $h1_words[1];
		my $h1_b1 = $h1_b_words[1];
		$is_paired_end = (($h11 =~ m/^1.+/) && (($h1_b1 =~ m/^2.+/))) ? 1 : 0;
	}
    }
    close(IN);
    return $is_paired_end;
}

# Separates a paired end fastq.  Assumes read 1 and read 2 are interleaved.
# 	@param $fq - paired end fastq file
# 	@param $read1_outfile - the output file name for read1
# 	@param $read2_outfile - the output file name for read2
# 	@returns - (nothing)
sub separatePairedEndFq {
    my ($fq, $read1_outfile, $read2_outfile) = @_;

    open(OUT_1, ">$read1_outfile");
    open(OUT_2, ">$read2_outfile");
    open(IN,    $fq) or die("Unable to open fq: $fq\n");
    while ( defined(my $h1_a = <IN>)
        and defined(my $s_a  = <IN>)
        and defined(my $h2_a = <IN>)
        and defined(my $q_a  = <IN>)
        and defined(my $h1_b = <IN>)
        and defined(my $s_b  = <IN>)
        and defined(my $h2_b = <IN>)
        and defined(my $q_b  = <IN>))
    {
        print OUT_1 $h1_a, $s_a, $h2_a, $q_a;
        print OUT_2 $h1_b, $s_b, $h2_b, $q_b;
    }
    close(IN);
    close(OUT_1);
    close(OUT_2);
}

# Converts a FASTQ file to FASTA and qual files
sub fastqToFastaQual {
    my ($fq, $fa, $qual) = @_;
    system("jgi_fastq_to_fasta_qual $fq $fa $qual");
}

##############################################
# 			         FASTA
##############################################

# Counts the nubmer of sequences in a given FASTA file
# 	@param $fa - the fasta file
# 	@returns - the number of sequences
sub cntSeqs {
    my ($fa) = @_;
    my $cnt = `grep -c '^>' $fa`;
    chomp($cnt);
    return $cnt;
}

# Selects the sequences in the given list of ids from a FASTA
# 	@param $ids_to_select_ref - an array ref of the ids to select
# 	@param $fasta - the input FASTA file to select from
# 	@param $outfile - the output file containing the selected sequences
# 	@returns - (nothing)
sub selectSeqsById($$$) {
    my ($ids_to_select_ref, $fasta, $outfile) = @_;

    # load ids into hash
    my %select_ids = ();
    foreach my $id (@$ids_to_select_ref) {
        $select_ids{$id} = 1;
    }

    # select
    my $in_ob  = Bio::SeqIO->new(-file => $fasta,      '-format' => 'Fasta');
    my $out_ob = Bio::SeqIO->new(-file => ">$outfile", -format   => 'fasta');
    my $num_selected = 0;
    while (my $seq = $in_ob->next_seq()) {
        if (defined($select_ids{ $seq->id })) {
            $out_ob->write_seq($seq);
            $num_selected++;
        }
    }
}

# Renames sequences with a simple integer, starting at 1
sub renameSeqs {
    my ($fasta, $outfile) = @_;
    my $in_ob  = Bio::SeqIO->new(-file => $fasta,      '-format' => 'Fasta');
    my $out_ob = Bio::SeqIO->new(-file => ">$outfile", -format   => 'fasta');
    my $cnt    = 1;
    while (my $seq = $in_ob->next_seq()) {
        my $seq_obj = Bio::Seq->new(
            -seq        => $seq->seq,
            -display_id => $cnt++,
            -desc       => $seq->desc,
            -alphabet   => "dna"
        );
        $out_ob->write_seq($seq_obj);
    }
}

# Splits a FASTA file
sub splitFa {
    my ($fasta, $num_out, $outfile_prefix) = @_;
    if ($num_out == 0) {
        die("Cannot split a file into 0 jobs.");
    }

    # determine how many reads we want in each file
    my $num_seqs = cntSeqs($fasta);
    my $num_seqs_per_file = int($num_seqs / $num_out) + ($num_seqs % $num_out == 0 ? 0 : 1);

    # split
    my $curr_file_cnt = 1;
    my $curr_seq_cnt  = 0;
    my $out_ob = Bio::SeqIO->new(-file => ">$outfile_prefix.$curr_file_cnt", -format => 'fasta');
    my @split_files = ("$outfile_prefix.$curr_file_cnt");
    my $in = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
    while (my $seq = $in->next_seq()) {
        if ($curr_seq_cnt >= $num_seqs_per_file) {
            $out_ob->close;
            $curr_file_cnt++;
            $curr_seq_cnt = 0;
            $out_ob =
              Bio::SeqIO->new(-file => ">$outfile_prefix.$curr_file_cnt", -format => 'fasta');
            push(@split_files, "$outfile_prefix.$curr_file_cnt");
        }
        
        $out_ob->write_seq($seq);
        $curr_seq_cnt++;
    }
    close(IN);
    $out_ob->close;

    if (@split_files != $num_out) {
        die(    "Meant to split into $num_out chunks, but really split into "
              . scalar @split_files
              . " chunks. Num per file: $num_seqs_per_file.");
    }

    return \@split_files;
}

sub splitFaRandomly {
    my ($fasta, $num, $outfile_prefix) = @_;
    if (!(-e $fasta)) {
        die("Input fasta file does not exist: $fasta");
    }
    if ($num > 1000) {
        print "WARNING: may have too many open file handles.\n";
    }

    # open output files
    my @split_files = ();
    my @out_fhs     = ();
    for (my $i = 1 ; $i <= $num ; $i++) {
        my $file_name = "$outfile_prefix.$i";
        my $out_ob = Bio::SeqIO->new(-file => ">$file_name", -format => 'fasta');
        push(@out_fhs,     $out_ob);
        push(@split_files, $file_name);
    }

    # randomly write to output files
    my $in = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
    while (my $seq = $in->next_seq()) {
        my $rand_file_index = int(rand($num));
        my $out_ob          = $out_fhs[$rand_file_index];
        $out_ob->write_seq($seq);
    }

    # make sure output files have something written
    foreach my $file (@split_files) {
        if (!(-e $file) || -s $file == 0) {
            die("Split file is empty or does not exist: $file\n");
        }
    }

    return (\@split_files);
}

# Randomly spliuts a FASTA file with sequences only occupying
# one line.
sub splitSingleLineFastaRandomly {
    my ($filename, $num_out, $pe) = @_;

    # open output files
    my @split_files = ();
    my @out_fhs     = ();
    for (my $file_cnt = 1 ; $file_cnt <= $num_out ; $file_cnt++) {
        my $split_file = "$filename.$file_cnt";
        push(@split_files, $split_file);
        my $fh;
        open($fh, ">$split_file") || die("Unable to open $split_file for writing\n");
        push(@out_fhs, $fh);
    }

    # write to several output files
    open(IN, $filename) || die("Unable to open input file: $filename\n");
    while (defined(my $header = <IN>) && defined(my $seq = <IN>)) {

        # randomly choose an output file
        my $out_fh = $out_fhs[ int(rand($num_out)) ];

        # write to the appropriate file (ensures every sequence is output somewhere)
        print $out_fh $header, $seq;
        if ($pe && defined(my $header2 = <IN>) && defined(my $seq2 = <IN>)) {
            print $out_fh $header2, $seq2;
        }
    }
    close(IN);

    # close output files
    for my $out_fh (@out_fhs) {
        close($out_fh);
    }

    return \@split_files;
}

# Separates a paired end fastq.  Assumes read 1 and read 2 are interleaved.
# 	@param $fq - paired end fastq file
# 	@param $read1_outfile - the output file name for read1
# 	@param $read2_outfile - the output file name for read2
# 	@returns - (nothing)
sub separatePairedEndFa {
    my ($fa, $read1_outfile, $read2_outfile) = @_;

    open(OUT_1, ">$read1_outfile")
      || die("Unable to open file for separation of reads: $read1_outfile\n");
    open(OUT_2, ">$read2_outfile")
      || die("Unable to open file for separation of reads: $read2_outfile\n");
    open(IN, $fa) or die("Unable to open fa: $fa\n");
    while ( defined(my $h1_a = <IN>)
        and defined(my $s_a  = <IN>)
        and defined(my $h1_b = <IN>)
        and defined(my $s_b  = <IN>))
    {
        print OUT_1 $h1_a, $s_a;
        print OUT_2 $h1_b, $s_b;
    }
    close(IN);
    close(OUT_1);
    close(OUT_2);
}

# Selects sequences from fasta with length >= min_len
# 	@param $in - the input fasta file
# 	@param $out - the output fasta file
# 	@param $min_len - the minimum length of a sequence to be included
# 	@returns - the number of contigs selected
sub filterByLen {
    my ($in, $out, $min_len) = @_;
    (-e $in) || die("Input file does not exist: $in");

    my $in_ob  = Bio::SeqIO->new(-file => $in,     '-format' => 'Fasta');
    my $out_ob = Bio::SeqIO->new(-file => ">$out", -format   => 'fasta');
    my $cnt    = 0;
    while (my $seq = $in_ob->next_seq()) {
        if (length($seq->seq()) >= $min_len) {
            my $seq_obj = Bio::Seq->new(
                -seq        => $seq->seq,
                -display_id => $seq->id,
                -desc       => $seq->desc,
                -alphabet   => "dna"
            );
            $out_ob->write_seq($seq_obj);
            $cnt++;
        }
    }
    return $cnt;
}

# Calculates common statistics for a given set of contigs.
# 	@param $fasta - the fasta file
# 	@returns - (num_seqs, total_len, n50, min, max)
sub calcContigStatsMinLen {
    my ($fasta, $min_len) = @_;

    # get sequence lengths
    my %seq_lengths = ();
    my $in_ob = Bio::SeqIO->new(-file => $fasta, '-format' => 'Fasta');
    while (my $seq = $in_ob->next_seq()) {
        $seq_lengths{ $seq->id } = length($seq->seq);
    }

    # length stats
    my @sorted_lengths_raw = sort { $a <=> $b } values %seq_lengths;
    my @sorted_lengths = ();
    my ($short_cnt, $med_cnt, $long_cnt) = (0, 0, 0);
    foreach my $length (@sorted_lengths_raw) {
        if ($length >= $min_len) {
            push(@sorted_lengths, $length);
            if ($length >= 1_000) {
                $long_cnt++;
            }
            elsif ($length >= 500) {
                $med_cnt++;
            }
            else {
                $short_cnt++;
            }
        }
    }

    # median
    my $median_len = median(\@sorted_lengths);

    my $total_length = 0;
    foreach my $len (@sorted_lengths) {
        $total_length += $len;
    }

    my $total_n50_len = 0;
    my $n50           = -1;
    for (my $i = $#sorted_lengths ; $i >= 0 ; $i--) {
        my $length = $sorted_lengths[$i];
        $total_n50_len += $length;
        if ($total_n50_len > ($total_length / 2)) {
            $n50 = $length;
            last;
        }
    }

    my $min      = $sorted_lengths[0];
    my $max      = $sorted_lengths[$#sorted_lengths];
    my $num_seqs = @sorted_lengths;

    return ($num_seqs, $total_length, $n50, $min, $max, $short_cnt, $med_cnt, $long_cnt,
        $median_len);
}

#################################################################################
# DEPRECATED (below) - DO NOT USE IN FUTURE CODE
#################################################################################

return 1;
