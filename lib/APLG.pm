package APLG;

=head1 NAME

APLG

=head1 SYNOPSIS

use lib "/home/klabutti/bin/modules";
use APLG;
my $run = shift @ARGV;
my $runObj = APLG->new("$run");

Ex usage:
#!/jgi/tools/bin/perl -w
use strict;
use lib "/home/klabutti/bin/modules";
use APLG;
 
my $run = shift @ARGV;

my $runObj = APLG->new("$run", "subdirName");

my $RUN = $runObj->getRun();

print "RUN: $RUN\n";


The available functions (pretty self explanatory):
getPre
getRef    
getData    
getRun    
getDate    
getVersion 


=head1 DESCRIPTION
This module parses an allpathsLG run for basic information about location, date, version.


=head1 AUTHOR(s)

Kurt M. LaButti

=head1 HISTORY
25oct2012 removed legacy  s/ASSEMBLY/ ASSEMBLY/ ..not sure what that was for..
31may2012 fixed parsing of assembly_stats file to incl lines w/o / @ end
28oct2011 added s///; functions to force a space in front of PRE, DATA, etc
28OCT2011 modified regexp to find REF and DATA
=over

=item *

Kurt M. LaButti 18october2011 Creation

=back

=cut

#============================================================================#
use strict;
use Carp;
use Carp qw(cluck);
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../lib";


my $_DEBUG = 0;
                                             
#============================================================================#
sub new {
    #assumes run/ASSEMBLIES/test is the dir structure    
    my $class       = shift;
    my $allPathsRun = shift; 
    my $subdir      = "test";
    my $fullPath    = "$allPathsRun/ASSEMBLIES/$subdir";
    my $statsFile   = "assembly_stats.report";


#test if run dir exists
    if ( !-d "$allPathsRun" ) {
        #confess "ERROR: APLG run $allPathsRun does not exist.\n";
	return "";
    }
#test if run dir subdirs exist
    if ( !-d "$fullPath" ) {
        confess "ERROR: APLG run subdir $fullPath does not exist.\n";
    }
#test if stats report file exists
   if ( !-s "$fullPath/$statsFile" ) {
        confess "ERROR: APLG stats file $fullPath/$statsFile does not exist.\n";
    }
    
    my $self = {};
    
    bless $self, $class;

    $self->_parseFile("$fullPath/$statsFile");
    
    return $self;

}

#============================================================================#
sub _parseFile {

    my $self = shift;    
    my $statsFile = shift;
    
    my ($assemDate, 
        $assemVersion,
        $PRE, 
        $REFERENCE_NAME, 
        $DATA,
        $RUN,
        $cmdLine);
    my %hash;
    # %hash = 
    #       {BLOCK_NAME} = value

           
    unless (open FH, $statsFile ) {
        confess "ERROR: failed to open file $statsFile: $!";
    }
 
    my $pidSeen = 0;    
    while (<FH>) {
        if ($_ =~ /pid\=/) {
	    ++$pidSeen;       
	    my @info = split /\s+/, $_;
            $assemDate = "$info[2]$info[1]$info[4]";
            $assemVersion = "$info[13]";
        } else {
	    last if ($_ =~ /^\-+/ && $pidSeen > 0);
            chomp;	    
            $_ =~ s/^\s+//;
            $_ =~ s/ \\//;
            $cmdLine .= $_;
        }
    }
    close FH;
        
    # note this does not catch 
    # the last bit of the command 
    # b/c the line does not have 
    # a "\" at the end.

    #make sure there is a space in front 
    $cmdLine =~ s/.*PRE/ PRE/;
    $cmdLine =~ s/RUN/ RUN/;
    $cmdLine =~ s/DATA/ DATA/;
    $cmdLine =~ s/RUN/ RUN/;
    $cmdLine =~ s/SUBDIR/ SUBDIR/;
    #$cmdLine =~ s/ASSEMBLY/ ASSEMBLY/;

    #print STDERR "CHK: [$cmdLine]\n";
    $cmdLine =~ /DATA=(\S+)/; $DATA = $1; 
    $cmdLine =~ /PRE=(\S+)/; $PRE = $1;
    $cmdLine =~ /RUN=(\S+)/; $RUN = $1;
    #print STDERR "DATA=$DATA\n";
    $DATA =~ /^((\w|\.|\_|\-)+)\/(\S+)/;
    #print STDERR "stuff=[$1] [$2] [$3]\n"; 
    $REFERENCE_NAME = $1;
    $DATA = $3;
    #print STDERR "D=[$DATA]\n";
    $assemDate =~ s/,//;
    #print STDERR "CHK2: [$cmdLine]\n";

    #print STDERR "CHECK\nPRE=$PRE\nDATA=[$DATA]\nREF=$REFERENCE_NAME\nRUN=$RUN\nDATE=$assemDate\nVER=$assemVersion\n";       
    $hash{'PRE'}      = $PRE;  
    $hash{'REF_NAME'} = $REFERENCE_NAME;
    $hash{'DATA'}     = $DATA;
    $hash{'RUN'}      = $RUN;
    $hash{'DATE'}     = $assemDate;    
    $hash{'VERS'}     = $assemVersion; 
    #print STDERR "CHECK\nPRE=$PRE\nDATA=$DATA\nREF=$REFERENCE_NAME\nRUN=$RUN\nDATE=$assemDate\nVER=$assemVersion\n";       

    $self->{_hash} = \%hash;
}



#============================================================================#
sub getPre {    
    my $self = shift;
    return $self->{_hash}{'PRE'};
}

sub getRef {    
    my $self = shift;
    return $self->{_hash}{'REF_NAME'};
}

sub getData {    
    my $self = shift;
    return $self->{_hash}{'DATA'};
}

sub getRun {    
    my $self = shift;
    return $self->{_hash}{'RUN'};
}

sub getDate {    
    my $self = shift;
    return $self->{_hash}{'DATE'};
}


sub getVersion {    
    my $self = shift;
    return $self->{_hash}{'VERS'};
}

1;


