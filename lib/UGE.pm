package UGE;

=head1 NAME

ktools UGE handler

=head1 SYNOPSIS


Ex usage:
#!/jgi/tools/bin/perl -w
use strict;
use lib "/global/u1/k/klabutti/lib";
use UGE;
 

my @jobs = ( .... )

my $obj = UGE->new(); #create obj
$obj->setSleep(60);  #set params
$obj->useUGE(@jobs); #run jobs


=head1 DESCRIPTION

=head1 AUTHOR(s)

Kurt M. LaButti

=head1 HISTORY

2012nov29 added -N jobName back in
2012oct11 added reset for _params in run sub

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
use File::Temp;
use File::Copy;
use File::Path;
use File::Path qw(rmtree mkpath);
use File::Basename;
use Cwd;
use Cwd qw (realpath getcwd abs_path);
use lib "/global/u1/k/klabutti/lib";

my $_DEBUG = 0;
                                             

#==============================================================================#
sub new {
    my $class = shift;
    my %params = @_;
    my $self = {};

#
#switches and check sleep time
    $self->{_noParams} = ''; #sets no params!
    $self->{_sleep}    = 1800;  #sleep 30 min in between checking jobs
#
#UGE parameter options, defaults set for typical fungal assem job
    $self->{_basic}    = '-m abe -cwd -V -R y -b yes -j yes -w e';
    $self->{_email}    = 'klabutti@lbl.gov';  
    $self->{_proj}     = 'fungal-assembly.p'; 
    $self->{_slots}    = '8';
    $self->{_ram}      = '5'; #in Gb
    $self->{_time}     = '10:00:00'; #in hr:min:sec
#    $self->{_jobName}  = 'no_name';
    $self->{_params}   = "$self->{_basic} "
	."-M $self->{_email} "
#	."-N $self->{_jobName} "
	."-P $self->{_proj} "
	."-l ram.c=$self->{_ram}G,h_rt=$self->{_time} "
	."-pe pe_slots $self->{_slots} ";
    
    bless $self, $class;
    return $self
}

#==============================================================================#
sub setBasicParams {
    my $self = shift;
    $self->{_basic} = shift;
}
sub setEmail {
    my $self = shift;
    $self->{_email} = shift;
}
sub setRam {
    my $self = shift;
    my $ram = shift;
    $self->{_ram} = $ram;
}
sub setSlots {
    my $self = shift;
    $self->{_slots} = shift;
}
sub setTime {
    my $self = shift;
    $self->{_time} = shift;
}
sub setJobName {
    my $self = shift;
    $self->{_jobName} = shift;
}
sub setProj {
    my $self = shift;
    $self->{_proj} = shift;
}
sub setNoParams {
    my $self = shift;
    $self->{_setNoParams} = 1;
}
 
sub setSleep {
   my $self = shift;
    $self->{_sleep} = shift;
}
#==============================================================================#
sub useUGE {
    my $self = shift;
    my @jobs = @_;
   


#ensure any updates to params are kept
    $self->{_params}   = "$self->{_basic} "
	."-M $self->{_email} "
	."-P $self->{_proj} "
	."-l ram.c=$self->{_ram}G,h_rt=$self->{_time} "
	."-pe pe_slots $self->{_slots} ";
    $self->{_params} .= " -N $self->{_jobName} " if ($self->{_jobName});
#    k::msg("params: $self->{_ram}");
    
    
    my %jobIds;
     
    foreach my $job (@jobs) {
        chomp $job;

        # construct qsub 
        my $qsubCmd = ($self->{_setNoParams})
	    ? "qsub $job" #must have  '  params "script"  '
	    : "qsub $self->{_params} \"$job\"";
	
        # send Qsub
	print STDERR "UGE: $qsubCmd\n";  
        my $jobInfo = `$qsubCmd`; 
        chomp $jobInfo;
        $jobInfo =~ /Your\s+job\s+(\d+)\s+\(\"(\S+)\"\)\s+has\s+been\s+submitted/;
        my $obsJobId   = $1;
        my $obsJobName = $2;
        ++$jobIds{$obsJobId};
        k::msg("jobID=$obsJobId jobName=$obsJobName cmd=$qsubCmd");    
    } 
     
#initial sleep to ensure jobs are in cached copy of qstat
    k::msg("sleeping for 240s to wait for cached qstat");
    sleep 240;

    #monitor UGE job id(s)
    monitorUGE($self,\%jobIds);
    k::msg("=all jobs have completed=");
} 
  
#==============================================================================#
sub monitorUGE {
    my $self = shift;
    my ($href) = @_;
    my %jobIds = %$href;
    my $continue = 0;
 
     
    while ($continue == 0) {
        my $numRunningJobs = 0;
	$continue = 1;
 
        foreach my $id (keys %jobIds) {
            chomp $id;

	    my @status = split /\s+/, `isjobcomplete $id`;


	    if ($#status == 1) {
		k::msg("$status[0] $status[1] still in the queue");
		++$numRunningJobs;
	    } else {
		k::msg("$status[0] not in queue, deleting");
		delete $jobIds{$id};
            }
        } 
 
	$continue = 0 if ($numRunningJobs > 0);
	
        if ($continue == 0) { 
            k::msg("$numRunningJobs jobs running ... sleep=$self->{_sleep}"); 
            sleep  $self->{_sleep};
        }
	
	
    } 
     
}  
 



1;


 
