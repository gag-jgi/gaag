package k;

=head1 NAME

ktools

=head1 SYNOPSIS

use lib "/global/u1/k/klabutti/lib";
use k;

Ex usage:
#!/jgi/tools/bin/perl -w
use strict;
use lib "/global/u1/k/klabutti/lib";
use k;
 

=head1 DESCRIPTION

=head1 AUTHOR(s)

Kurt M. LaButti

=head1 HISTORY

#16oct2012 added libinfo function
#10oct2012 updated assemStats regexp

=over

=item *
Kurt M. LaButti 14march2012 Creation

=back

=cut

#============================================================================#

use Carp;
use Carp qw(cluck);
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../lib";
use Bio::SeqIO;
use Bio::Seq;
use File::Temp;
use File::Copy;
use File::Path;
use File::Path qw(rmtree mkpath);
use File::Basename;
use Cwd;
use Cwd qw (realpath getcwd abs_path);
use lib "/global/u1/k/klabutti/lib";
use lib "/house/groupdirs/QAQC/scripts/versions/DEFAULT/lib";
use PGF::Utilities::JobManager;
use List::Util qw(sum);

my $_DEBUG = 0;




#==============================================================================#
sub gcsw { 
    my ($fasta,$shredLen,$ratchet) = @_;
    
    my $overlap = $shredLen - $ratchet; 
    my $head;
#    k::msg("sub_shred: $header shredLen=$shredLen shredOl=$overlap fa=$fasta") if ($optDebug);
    my $print;
    my $beg  = 0;
    my $x    = 0;
    my $stop = 0;
    my $checkBeg = $beg;
    my $checkEnd = $beg+$shredLen;
    my $checkShredLen = $checkEnd - $checkBeg;
    my $seqLen = length($fasta);
    my @gcs;
#
# main block
    if ($seqLen >= $shredLen) {
        #seqLen is > shredLen, shred er up
        while ($checkBeg < $seqLen && $stop == 0) {# && $checkEnd <= $seqLen) {
            #beg,end are within seq bounds
            if ($checkBeg < $seqLen && $checkEnd < $seqLen) {
		#       print STDERR "##SUBSTR [within] $head b=$checkBeg e=$checkEnd l=$seqLen substr($checkBeg ,$shredLen);\n" if ($optDebug);
                $print = substr($fasta, $checkBeg ,$shredLen);
                #end is past seq bounds, but stop if the tot shred size is short
            } elsif ($checkBeg < $seqLen && $checkEnd >= $seqLen) {
		#       print STDERR "##SUBSTR [end past] $head $checkBeg $checkEnd [beg > $seqLen]\n" if ($optDebug);
                $print =  substr($fasta, $checkBeg);
                ++$stop;
            } else {
                #no more seq, do nothing
                print STDERR "SKIP [no more seq]\n" if ($optDebug);
                last;        
            }
	    
            my $gc = fullgc($print);
            push(@gcs,$gc);
            chomp $print if ($print);
	    #           $print =~  s/(.{50})/$1\n/g if ($seqLen > 50 && $print);
#           print STDOUT ">${header}_${x}_$gc\n$print\n" if ($print); #beg=$checkBeg end=$checkEnd origLen=$seqLen shredLen=$checkShredLen ratchet=$overlap (ol=$ratchet)
            $checkBeg += $overlap;
            $checkEnd += $overlap;
            ++$x;        
        }
	
    } else {
	
        chomp $fasta if ($fasta);
        my $gc = fullgc($fasta);
        push(@gcs,$gc);
	
#       print STDOUT ">${header}_0_$gc\n$fasta\n" if ($fasta);
    }
    
#   print join "\n", @gcs;
#   print "\n";
    my $agc = sprintf "%.2f", mean(@gcs);
    return $agc;
}

#==============================================================================#
sub mean {
    return sum(@_)/@_;
}

#==============================================================================#
sub fullgc {  
    my ($seq) = @_;
    my $gc = ($seq =~ tr/[GgCc]//);
    my $ln = length($seq);
    my $result = (($gc/$ln)*100);
}  
       
#==============================================================================#
sub gc { 
    my ($seq) = @_;
    my $gc = ($seq =~ tr/[GgCc]//);
    my $ln = length($seq);
    my $result = sprintf "%.1f" ,(($gc/$ln)*100);
} 
                                             
#==============================================================================#
sub msg {
    my $msg = shift @_;
    my $date = setDate();
    print STDERR "[$date] $msg ...\n";
}

#==============================================================================#
sub setDateOld {
    my $date = `date`;
    chomp $date;
    return $date;
}
 
#==============================================================================#
sub setDate {
#
# Returns current date and time in format: 'MM/DD/YYYY HH24:MI:SS'
    my ($sec,$min,$hour,$day,$mon,$year) = localtime(time); 

#
# $year contains no. of years since 1900, to add 1900 to make Y2K compliant
    $year+=1900;  
#
# month abbreviations #$mon++; ZERO BASED, handling this w/ @abbr below
    my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );

    my $dmy  = sprintf( "%02d/%s/%04d", $day, $abbr[$mon], $year);
    my $time = sprintf( "%02d:%02d:%02d", $hour, $min, $sec);

    return "$dmy $time"; #example "14/Mar/2012 15:34:28"
}

#==============================================================================#
sub argv2list { #@ARGV input, distills paths, lists of paths. returns paths
    my $a_ref = shift @_;
    my @list;
    
    foreach my $input (@$a_ref) {
        $input = realpath($input);
        chomp $input;    
	k::msg("INPUT $input");
        #
        #is input a directory
        if (-d "$input") {
            # input is dir, push to array
	    k::msg("DIRECTORY");
            push(@list, "$input")
                if (-e "$input");
        } elsif (-f "$input") {
	k::msg("FILE/FOF");
            # input is a file (hopefully a list of paths)
            open FOF, "$input" or die "Can't open fof [$input] $!\n";
            while (my $el = <FOF>) {
                chomp $el;
                $el = realpath($el);
                #push dir onto array if it exists
                push(@list, "$el") if (-e "$el");
            } 
            close FOF;      
        } else { 
            usage("ERROR: Issue with input: $input");
        } 
	
    }  
    # return list of paths
    return \@list;
}

#==============================================================================#
sub check {
    my ($thing) = shift @_;
    my $result = 0;

#check for file>0 and dir
    if (-s $thing) {++$result;}
    if (-d $thing) {++$result;}

    if ($result == 0) {
	confess "FAIL! file or dir [$thing] does not exist or is size 0";
    } else {
	print STDERR "ok! [$thing]\n";
    }
}

#==============================================================================#
sub assemStats {
    my $assemFa = shift @_;
    my $num     = shift @_;
    my $debugPrint = shift @_;
    my $exe  = "/house/homedirs/j/jschmutz/fasta_stats2.linux";
    my $numN = ($num && $num =~ /\d+/)
	? $num
	: "1";
    confess "you noob! file [$assemFa] doesn't exist!" 
	unless (-s $assemFa);
    confess "you noob! file [$exe] doesn't exist!" 
	unless (-s $exe);


    print STDERR "NUMN to split on = $numN\n";


    my $output = `$exe -n $numN $assemFa`;
    my @tmp    = split /\n/, $output;


    if ($debugPrint) {
    for (0..$#tmp) {
	print "$_ $tmp[$_]\n";
    }
}

    $tmp[0] =~ /\:\s+(\d+)/;
    my $numS = $1;
    $tmp[1] =~ /\:\s+(\d+)/;
    my $numC = $1;
    $tmp[2] =~ /\:\s+(\S+\s+(\w+)?)/;
    my $ttlS = $1;
    $tmp[3] =~ /\:\s+(\S+\s+(\w+)?)/;
    my $ttlC = $1;
    $tmp[4] =~ /\:\s+(\d+\/\d+\S+(\s+\w+)?)/;
    my $n50S = $1;
    $tmp[5] =~ /\:\s+(\d+\/\d+\S+(\s+\w+)?)/;
    my $n50C = $1;
    
    my $head = "#numScaff,numCtg,ttlScaff,ttlCtg,scaffN50,ctgN50";
    return "$numS,$numC,$ttlS,$ttlC,$n50S,$n50C";
    

}

#==============================================================================#
sub guessFastqFormat {
    my ($fastq) = shift @_;
    chomp $fastq;
    confess "ERROR! file [$fastq] doesn't exist!" unless (-s $fastq);
    my $format = `/global/u1/k/klabutti/bin/fastq_guess_format $fastq`;
    chomp $format;
    return $format;
}

#==============================================================================#
sub gimmedir {
    my ($dir) = @_;

my @fullPathFiles;
   
    opendir DH, $dir or die "Cannot open $somedir: $!";
    my @files = grep { ! -d } readdir DH;
    closedir DH;

    foreach my $f (@files) {
	chomp $f;
	push(@fullPathFiles, "$dir/$f");
    }

    return \@fullPathFiles;
}

#==============================================================================#
sub gimmefulldir {
    my ($dir) = @_;
 
my @fullPathFiles;
   
    opendir DH, $dir or die "Cannot open $somedir: $!";
    my @files = readdir DH;
    closedir DH;
 
    foreach my $f (@files) {
        chomp $f;
        push(@fullPathFiles, "$dir/$f");
    }
 
    return \@fullPathFiles;
}

#==============================================================================#
sub getMd5 {
    my $thing = shift @_;
    my $thingToTag = realpath($thing);
    my $exe       =  "/jgi/tools/bin/md5sum-lite";
    unless (-s $thing) {
	confess "thing [$thing] does not exist or is empty!\n";
    }
    unless (-X $exe) {
	confess "exe [$exe] is not executeable by you!\n";
    }
    chomp $thing;   
    my $info = `$exe $thing`; 
    chomp $info;
    my ($tag,$loc) = split /\s+/, $info;

    return $tag;
}

#==============================================================================#
sub makeMd5 {
    my $el  = shift @_;
    my $exe =  "/jgi/tools/bin/md5sum-lite";
    my @out;

    my $fullPath = realpath($el);
    unless (-s $fullPath) {
	warn "[$fullPath] is size zero, skipping";
    }
    my ($name, $loc) = fileparse($fullPath);
#    $name =~ /(.*)\.\S+$/;
 #   $out = $1;
    $out = "${fullPath}.md5";
    
    my $cmd = "$exe $fullPath > $out";

    system($cmd);
   
    return $out if (-s $out);
}


1;
