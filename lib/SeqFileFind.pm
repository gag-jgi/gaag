package SeqFileFind;

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
19nov2012 added exit for null goiLoc path 
16oct2012 added exclusion for project type RNA, or library type pe_str
20sept2012 added ignore for "Transcriptome" in the organism name checker
14sept2012 added 	$split[0] =~ s/ $//; b/c of inconsistent entries in db
17aug2012 seq_file_find was modified to pull from ITS instead of RQC, so had to modify this module to sync the header key names
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
use Cwd qw (realpath getcwd abs_path);

my $_DEBUG = 0;
         
#warn "SFF is successfully loaded!\n";

                                    
#============================================================================#
sub new {
    #assumes run/ASSEMBLIES/test is the dir structure    
    my $class    = shift;
    my $organism = shift;

    #check org name
    if ( $organism !~  /\w+/ ) {
        confess "ERROR: organism must be a word\n";
    }
    
    my $self = {};
    
    bless $self, $class;

    $self->_getSffInfo("$organism");
    
    return $self;

}

#============================================================================#
sub _getSffInfo {
    my $self     = shift;    
    my $organism = shift;
    my %info;
#
#$info{$c}{all}        = "full,data,line"
#         {headerName} = value
#
    $self->_bio2data("$organism");
    
    return $self;
}


    
#============================================================================#
sub _bio2data {
    my $self         = shift;
    my $organism     = shift;
    my $SEQFILEFIND  = "/house/groupdirs/QAQC/scripts/seq_file_find";
    my %numOrganisms;  #$numOrganisms{org_name} = num

    unless (-x $SEQFILEFIND) {
	confess "ERROR: seq_file_find is not exec!\n"; 
    }

    my @seqData     = `$SEQFILEFIND -bio \"$organism\"`;

    unless (@seqData > 0) {
	confess "ERROR: no data found for $organism!\n"; 
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
# change isnert "kbp" to "kb"
	$line =~ s/kbp,/kb,/;
	#
	# store full line
	$info{$c}{'all'} = $line;
	#
	# split line
	my @split  = split /\,/, $line; 
	$split[0] =~ s/ $//;
	#
	# track number of unique org names seen
	++$numOrganisms{$split[0]} unless ($split[0] =~ /Transcriptome/i);
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
	confess "No data found for organism $organism: $!\n";
    }
    #
    # fail if >1 organism found
    if ((scalar keys %numOrganisms) > 1) {
	foreach my $o (keys %numOrganisms) {warn "organism=[$o]\n";}
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
    my $genomicEnd      = 0;
    my $organismName    = $hoh{0}{'ncbi_name'}; #was bioclass prior to ITS
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

	print "---->  $hoh{$num}{file_path}\n";
	#
	# skip cDNA lib types
	next if (exists $hoh{$num}{product_type} 
		 && $hoh{$num}{product_type} =~ /annotation/i);
	next if ($hoh{$num}{seq_project_name} =~ /RNA/ 
		 || $hoh{$num}{library_type} =~ /pe-str/);
	
	#
	# set variables 
	my $filePath = ($hoh{$num}{file_path});
	my $readLen  = ($hoh{$num}{run_profile} ne "na")
	    ? "$hoh{$num}{run_profile}"
	    : "0"; 
	my $fullInsert = ($hoh{$num}{lab_insert_size} =~ /\d+/) 
	    ? "$hoh{$num}{lab_insert_size}"
	    : "0";
	my $insert = $fullInsert;

#	print STDERR "FULL-INSERT: $fillInsert\n";
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
	my $name    = $hoh{$num}{library}; 
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



#============================================================================#
sub getQCCsv {
    my $self = shift;
    my $hoh_ref = $self{_hoh};
    my %hoh     = %$hoh_ref;

    my @data;
    my %tmp;
    my $libs;
    my $groups;
    my $readOri;
    my $genomicEnd      = 0;
    my $organismName    = $hoh{0}{'ncbi_name'}; #was bioclass prior to ITS
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
	next if (exists $hoh{$num}{product_type} 
		 && $hoh{$num}{product_type} =~ /annotation/i);
	next if ($hoh{$num}{seq_project_name} =~ /RNA/ 
		 || $hoh{$num}{library_type} =~ /pe-str/);


	#
	# lib function -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
	my ($readmeFullPath,
	    $estGenomeSize,
	    $insert, 
	    $stdDev, 
	    $mitoLoc, 
	    $goiLoc) = libinfo($hoh{$num}{library}); 
        #-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

	confess "goiLoc issue!! [$goiLoc]" unless ($goiLoc =~ /S+/);

	#
	# set variables 
	my $filePath = ($hoh{$num}{file_path});
	my $readLen  = ($hoh{$num}{run_profile} ne "na")
	    ? "$hoh{$num}{run_profile}"
	    : "0"; 

	#
	# set type according to insert
	 #print "INSERT! [$insert]\n";
	my $type = ($insert == 0)  
	    ? "long"  #pacbio
	    : ($insert < 1000)  
	    ? "fragment"  
	    : ($insert >= 1000)
	    ? "jumping" 
	    : "unknown";
	
	#
	# handle jump orientation
	my $jumpOri;
	if ($type eq "jumping") {
	    $jumpOri = ($hoh{$num}{library_type} =~ /LFPE/i)
		? "outward"
		: "inward";
	    
	    $qcOri = ($mitoLoc =~ /lfpe\_outward/) 
		? "outward"
		: "inward";
	
	    print STDERR "jump=$jumpOri [$hoh{$num}{library_type}] qc=$qcOri [$readmeFullPath]\n";
    
	    unless ($qcOri eq $jumpOri) { confess "JUMP ORIENTATION ISSUE\n";}
	}
	
	#
	# set libstats according to type
	my $libStats = ($type =~ /fragment/) 
	    ? "1,$insert,$stdDev,,,inward,0,0"
	    : ($type =~ /jumping/)
	    ? "1,,,$insert,$stdDev,$jumpOri,0,$genomicEnd"
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
	my $name    = $hoh{$num}{library}; 

	my $libName = ($name =~ /^PB\d+/)
	    ? "pacbio"
	    : "${name}_${insert}_${readLen}";
	# print "libNAME: [$libName]\n";	
	# IHPA_4kb_2x100_2155.2.1777
	my $groupName = "$libName" . "_$rundetails";
	
#
# print groups info
	$groups   .= sprintf "%-30s\t,%-20s\t,%s\n", 
	$groupName,$libName,$goiLoc;
	$groupName = "$libName" . "_mito" unless ($mitoLoc =~ /^na$/);
	$groups   .= sprintf "%-30s\t,%-20s\t,%s\n", 
	($groupName,$libName,$mitoLoc) unless ($mitoLoc =~ /^na$/);
	
	$tmp{$libName} = "${type}_${libStats}";

    }
#
# print libs info
    foreach my $libName (keys %tmp) {
#	print "LIBNAME=$libName=$tmp{$libName}\n";
	my ($lib,$readLen,$insert) = split /\_/, $libName;
	my ($type,$libStats)       = split /\_/, $tmp{$libName};
	$libs .= sprintf "%-20s\t,%s,%s,%s,%s\n", 
	$libName,"aplg",$orgNameNoSpaces,$type,$libStats;	
    }
    
    
    return "$groups","$libs"
}

#===============================================================#
sub libinfo {
#    my $self = shift;
    my $lib   = shift;
    my $qcdir = "/house/groupdirs/fungal_analysis/projects/ldl/";
    
    
#   print STDERR "lib=[$lib]\n";
 
   if ($lib =~ /\w\w\w\w/ && -d "$qcdir/$lib") {
       chdir "$qcdir/$lib";
       my $readmeLoc = `find -name "README.${lib}"`;
       chomp $readmeLoc;
 #      print STDERR "readme=[$readmeLoc]\n";
       my $readmeFullPath = realpath($readmeLoc);
       chomp $readmeFullPath;
 #      print STDERR "readmeFull=[$readmeFullPath]\n";
       my $stdDev  = "na";
       my $insert  = "na";
       my $estGenomeSize = "na";
       my $mitoloc = "na";
       my $goiLoc  = "na";
 
       #insert
       my $infoLine = 
           `grep weighted $readmeFullPath`; 
       my @tmp = split /\s+/, $infoLine;
       $stdDev = sprintf "%.0f", $tmp[-1];
       $insert = sprintf "%.2f", $tmp[-3];
 
#       print STDERR "line=[$infoLine]\ninsert=[$insert] stddev=[$stdDev]\n";
 
       #postQC file locs
       $mitoLoc =  
           `grep fastq $readmeFullPath | tail -2 | grep /mito`; 
       $goiLoc =  
           `grep fastq $readmeFullPath | tail -2 | grep -v /mito`; 
              
       #est genome size
       my $estSizeLine = 
           `grep "Main genome scaffold sequence total:" $readmeFullPath`; 
       @tmp = split /\s+/, $estSizeLine;

       if ($estSizeLine) {
	   
	   if ($#tmp == 6) {
	       my $multiplier = $tmp[6];
	       my $abbrevSize = $tmp[5];
	       $estGenomeSize = ($multiplier =~ /MB/)
		   ? ($abbrevSize * 1000000)
		   : ($multiplier =~ /KB/)
		   ? ($abbrevSize * 1000)
		   : 'idk';
	       
	   } elsif  ($#tmp == 5) {
	       $estGenomeSize = $tmp[5];
	   }
       }
       
       chomp($readmeFullPath,$estGenomeSize,$insert,$stdDev,$mitoLoc,$goiLoc);
       
       print "return $readmeFullPath,$estGenomeSize,$insert, $stdDev, $mitoLoc, $goiLoc\n"; 
      return $readmeFullPath,$estGenomeSize,$insert, $stdDev, $mitoLoc, $goiLoc;
 
   } else { 
 
       return "404";
 
   } 
    
    
} 



1;
