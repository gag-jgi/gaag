package SeqFileFindLib;

=head1 NAME

SeqFileFind  "seq file find"

=head1 SYNOPSIS

use lib "/house/homedirs/k/klabutti/lib";
use SeqFileFind;
my $bio    = shift @ARGV;
my $pmoObj = SeqFileFind->new("$bio");

Ex usage:
#!/jgi/tools/bin/perl -w
use strict;
use lib "/home/klabutti/bin/modules";
use SeqFilefind;
my $bioName = shift @ARGV;
my $sffObj  = SFF->new("$bioName"); 

my $array_ref = $sffObj->getAllData();
my @tmp = @$array_ref;
print "getAllData:\n";
print join "\n", @tmp;
print "\n";
 
print "getCsv:\n";
my ($group,$lib) = $sffObj->getCsv();
print "$group\n";
print "$lib\n";

 

The available functions (pretty self explanatory):
getAlldata
getCsv

=head1 DESCRIPTION



=head1 AUTHOR(s)

Kurt M. LaButti

=head1 HISTORY

03april2012 creation

=over

=item *

Kurt M. LaButti 

=back

=cut

#============================================================================#
use DBI;
use Carp;
use Carp qw(cluck);
use File::Basename;
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../lib";


my $_DEBUG = 0;
         
warn "SFF is successfully loaded!\n";

                                    
#============================================================================#
sub new {
    #assumes run/ASSEMBLIES/test is the dir structure    
    my $class    = shift;
    my $lib = shift;

    #check org name
    if ( $lib !~  /\w+/ ) {
        confess "ERROR: organism must be a word\n";
    }
    
    my $self = {};
    
    bless $self, $class;

    $self->_getSffInfo("$lib");
    
    return $self;

}

#============================================================================#
sub _getSffInfo {
    my $self     = shift;    
    my $lib = shift;
    my %info;
#
#$info{$c}{all}        = "full,data,line"
#         {headerName} = value
#
    $self->_bio2data("$lib");
    
    return $self;
}


    
#============================================================================#
sub _bio2data {
    my $self         = shift;
    my $lib     = shift;
    my $SEQFILEFIND  = "/house/groupdirs/QAQC/scripts/seq_file_find";
    my %numOrganisms;  #$numOrganisms{org_name} = num

    unless (-x $SEQFILEFIND) {
	confess "ERROR: seq_file_find is not exec!\n"; 
    }

    my @seqData     = `$SEQFILEFIND -lib \"$lib\"`;

    unless (@seqData > 0) {
	confess "ERROR: no data found for $lib!\n"; 
    }

    #
    # data -> hashOfhash
    my @header;
    my $c = 0;
    foreach my $line (@seqData) {
	chomp $line;
	#
	# grab header
	if ($line =~ /^#/) {
	    my @rawHeader = split /\,/, $line;
	    foreach my $head (@rawHeader) {
		$head =~ s/^\#//;
		$head =~ s/\d+\=//;
		push(@header, $head);
	    }
	    next;
	}
	#
	# store full line
	$info{$c}{'all'} = $line;
	#
	# split line
	my @split  = split /\,/, $line; 
	#
	# track number of unique org names seen
	++$numOrganisms{$split[0]};
	#
	# store each value under the header name
	for (0..$#header) {
	    $info{$c}{$header[$_]} = ($split[$_] ne "")
		? $split[$_]
		: "na";	
	}
	++$c;
    }	
    #
    # fail if no data
    my $numKeys = scalar(keys %info);
    if ($numKeys < 1) {
	confess "No data found for organism $lib: $!\n";
    }
    #
    # fail if >1 organism found
    if ((scalar keys %numOrganisms) > 1) {
	foreach my $o (keys %numOrganisms) {warn "organism=$o\n";}
	confess "more than 1 organism found in output, exiting: $!\n";
    }
    
    
    $self{_hoh} = \%info;
}

#============================================================================#
sub getAllData {
    my $self = shift;
    my @data;
    my $hoh_ref = $self{_hoh};
    my %hoh = %$hoh_ref;

    foreach my $num (keys %hoh) {
	push(@data, "$hoh{$num}{'all'}");
    }

    my @sort = sort @data;
    return \@sort;
}

#============================================================================#
sub getCsv {
    my $self = shift;
    my $hoh_ref = $self{_hoh};
    my %hoh     = %$hoh_ref;

    my @data;
    my %tmp;
    my $libs;
    my $groups;
    my $genomicEnd      = 44;
    my $organismName    = $hoh{0}{'bioclass'};
    my $orgNameNoSpaces = $organismName;
    $orgNameNoSpaces    =~ s/\s+/_/g;
    
    $libs = "library_name, "
	."project_name, "
	."organism_name, "
	."type,paired, "
	."frag_size, "
	."frag_stddev, "
	."insert_size, "
	."insert_stddev, "
	."read_orientation, "
	."genomic_start, "
	."genomic_end\n";
    
    $groups = "group_name, library_name, file_name\n";
    
    foreach my $num (sort {$a <=> $b} keys %hoh) {
	#
	# skip cDNA lib types
	next if (exists $hoh{$num}{lib_name} 
		 && $hoh{$num}{lib_name} =~ /^C/);
	#
	# set variables 
	my $filePath = ($hoh{$num}{file_path});
	my $readLen  = ($hoh{$num}{read_length} ne "na")
	    ? "$hoh{$num}{read_length}"
	    : "0"; 
	my $fullInsert = ($hoh{$num}{expected_insert_size} =~ /\d+/) 
	    ? "$hoh{$num}{expected_insert_size}"
	    : "0";
	my $insert = $fullInsert;
	$insert    =~ s/kb//i;
	$insert    = ($insert * 1000);
	#
	# set type according to insert
	my $type = ($insert == 0)  
	    ? "long"  #pacbio
	    : ($insert < 1000)  
	    ? "fragment"  
	    : ($insert >= 1000)
	    ? "jumping" 
	    : "unknown";
	#
	# set stdev based on type	
	my $stdDev   = ($type =~ /fragment/) 
	               ? 50 
		       : ($type =~ /jumping/) 
	               ? 500
		       : '';
	#
	# set libstats according to type
	my $libStats = ($type =~ /fragment/) 
	    ? "1,$insert,$stdDev,,,inward,0,0"
	    : ($type =~ /jumping/)
	    ? "1,,,$insert,$stdDev,inward,0,$genomicEnd"
	    : ($type =~ /long/)
	    ? "0,,,,,,," 
	    : "unknown_lib";
	#
	# filename to details
	# example: 2265.6.1840.GGCTAC.fastq.gz
	# pacbio example: filtered_subreads.fastq
	my ($filename, $path) = fileparse("$filePath");
	my $rundetails;
	chomp $filename;
	
	if ($path =~ /pacbio\/jobs\/\d+\/(\d+)\/data/) {
	    #if pacbio stick on the number before data 	
	    $rundetails = $1;
	} else {
	    #if not pacbio then use the flowcell info
	    $filename =~ /^((\d+\.){3})(\S+\.)?fastq(\.gz)?$/;
	    $rundetails = $1;
	    $rundetails =~ s/\.$//;
	}
	#
	# create libname
	# example: IHPA_4kb_2x100 or pacbio_018595
	my $name    = $hoh{$num}{lib_name}; 
	my $libName = ($name =~ /^PB\d+/)
	    ? "pacbio"
	    : "${name}_${fullInsert}_${readLen}";
	
	# IHPA_4kb_2x100_2155.2.1777
	my $groupName = "$libName" . "_$rundetails";
	
	$groups .= sprintf "%-30s\t,%-20s\t,%s\n", $groupName,$libName,$filePath;
	
        $tmp{$libName} = "${type}_${libStats}";
    }
    foreach my $libName (keys %tmp) {
	my ($lib,$readLen,$insert) = split /\_/, $libName;
	my ($type,$libStats)       = split /\_/, $tmp{$libName};
	$libs .= sprintf "%-20s\t,%s,%s,%s,%s\n", $libName,"aplg",$orgNameNoSpaces,$type,$libStats;
    }    
    return "$groups","$libs"
}


1;


