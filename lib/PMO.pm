package PMO;

=head1 NAME

PMO

=head1 SYNOPSIS

use lib "/global/u1/k/klabutti/lib";
use PMO;
my $projectId= shift @ARGV;
my $pmoObj = PMO->new("$projectId");

Ex usage:
#!/jgi/tools/bin/perl -w
use strict;
use lib "/home/klabutti/bin/modules";
use PMO;
 
my $projId = shift @ARGV;

my $pmoObj = APLG->new("$projId");

my $piName = $pmoObj->getPiName();

print "PI name: $PiName\n";


The available functions (pretty self explanatory):
getFinalSubmissionEnd   
getSPID   
getNcbiOrgName   
getPiEmail   
getSpecies   
getDraftSubmissionEnd   
getEmbargoEndDate   
getDraftSubmissionStart
getNcbiTaxId   
getPiName   
getProgramYear   
getPmoProjId   
getFinalSumbissionStart   
getGenus
getPiHtml   
getStrain   
getOrganismName   
getOrganismNameNoSpaces   
getDataReleaseControl  

=head1 DESCRIPTION
This module accesses pmo database for project assoc. info given a project id


=head1 AUTHOR(s)

Kurt M. LaButti, Brian Foster (pmo.pl that was modifed to make this module)

=head1 HISTORY

01march2012 added organismNameNoSpaces
21dec2011  added confess if $bio doesn't exist (meaning projId returned zilch)

=over

=item *

Kurt M. LaButti 15December2011 Creation

=back

=cut

#============================================================================#
use DBI;
use Carp;
use Carp qw(cluck);
use FindBin qw($RealBin);
use lib "$FindBin::RealBin/../lib";


my $_DEBUG = 0;
                                             
#============================================================================#
sub new {
    #assumes run/ASSEMBLIES/test is the dir structure    
    my $class       = shift;
    my $projectId   = shift; 

    #test if projId is 7 digi
    if ( $projectId !~  /\d{7}/ ) {
        confess "ERROR: projectId must be a 7 digit number\n";
    }
    
    my $self = {};
    
    bless $self, $class;

    $self->_getPmoDbInfo("$projectId");
    
    return $self;

}

#============================================================================#
sub _getPmoDbInfo {

    my $self = shift;    
    my $projectId = shift;
    my %info;
    my %hash_ref;
 
    $ENV{'ORACLE_HOME'} = "/jgi/tools/oracle_client/DEFAULT";
    $ENV{'TNS_ADMIN'}   = $ENV{ORACLE_HOME} . "/bin";


    my ($cursor,$dbh, $sql);
    $dbh = DBI->connect('dbi:Oracle:psf.world','siagro','i48ryd6');
    
    $sql = q
    {
        SELECT
	    s.directory_number project_number,
            s.bio_classification_name
        FROM
                projects.dt_gpts_dna_samples s
	WHERE s.directory_number  =  ? 
    };

    $cursor = $dbh->prepare("$sql");
    $cursor->execute($projectId);

    my $hash_ref = $cursor->fetchall_hashref("PROJECT_NUMBER");
    use Data::Dumper;
    $cursor->finish;
    my $bio = $hash_ref->{$projectId}{BIO_CLASSIFICATION_NAME};

    #exit if no bio
    unless($bio) {confess "projId does not exist!\n";}
   
    $sql = q
    {
	select proposals.contact_api.get_contact_name(pi_contact_id,'last', 'Y', 'N') "PI",
	pmo_project_id "PMO Project ID",
	PMO_PROJECT_NAME "Name",
	program_and_year "Program/Year",
	embargo_end_date as "Embargo End Date",
	data_release_control as "Data Release Control",
	ncbi_organism_name "NCBI Organism Name",
	ncbi_tax_id "NCBI Tax ID",
	draft_asm_gbk_start_date "Draft Submission Start",
	draft_asm_gbk_end_date "Draft Submission End",
	final_asm_gbk_start_date "Final Submission Start",
	final_asm_gbk_end_date "Final Submission End",
	genus "Genus",
	species "Species",
	strain "Strain",
	spid "SPID"
	    from proposals.dt_pmo_projects_report p
    }; 
    $cursor = $dbh->prepare( "$sql");
    $cursor->execute();
    my $array_ref = $cursor->fetchall_arrayref();
    


# 
#cleanup names
# 
    my %name2info;
    my $name; 
    my $record;

    for ($i=0;$i<@{$array_ref};$i++){
	for  (my $j=0;$j<@{$array_ref->[$i]};$j++){
	    if (! defined ${$array_ref->[$i]}[$j]){
		${$array_ref->[$i]}[$j] = "";
	    }
	}
	my %temphash = ();
	($temphash{"PI"},
	 $temphash{"PMO_PROJECT_ID"},
	 $temphash{"Name"},
	 $temphash{"Program/Year"},
	 $temphash{"Embargo_End_Date"},
	 $temphash{"Data_Release_Control"},
	 $temphash{"NCBI_Organism_Name"},
	 $temphash{"NCBI_Tax_ID"},
	 $temphash{"Draft_Submission_Start"},
	 $temphash{"Draft_Submission_End"},
	 $temphash{"Final_Submission_Start"},
	 $temphash{"Final_Submission_End"},
	 $temphash{"Genus"},
	 $temphash{"Species"},
	 $temphash{"Strain"},
	 $temphash{"SPID"}) = @{$array_ref->[$i]};

	if (! $temphash{'Name'} || $temphash{"Name"} =~ /^Placeholder/){
	    next;
	}
	if ($temphash{"PI"}=~/^.*mailto\:(.*?)\">(.*)\<.*$/){
	    $temphash{"PI_NAME"}=$2;
	    $temphash{"PI_EMAIL"}=$1;
	}
	else{
#	warn "no mail for $temphash{PI}\n";
	    next;
	}
	
	$temphash{'Name'} =~ s/\s+/ /g;
	$temphash{'Name'} =~ s/\s+$//;
	$temphash{'Name'} =~ s/^\s+//;
	push @{$name2info{$temphash{"Name"}}},\%temphash;
	
    }
    
if (scalar(@{$name2info{$bio}}) > 1){
    warn "XXXXXXX multiple entries in pmo database for $bio XXXXXXX\nCheck by hand\nOnly reporting first contact info\n";
}
    foreach my $key (keys %{$name2info{$bio}[0]}){
#	print "$key=$name2info{$bio}[0]->{$key}\n";
	$info{$key} = "$name2info{$bio}[0]->{$key}";
    }
    
        
    $cursor->finish();
    $dbh->disconnect();
    
    
    $self->{_hash} = \%info;
}



#============================================================================#
sub getFinalSubmissionEnd {    
    my $self = shift;
    
    return  ( exists $self->{_hash}{Final_Submission_End} ) ?
	$self->{_hash}{Final_Submission_End} :
	"na";
}
sub getSPID {    
    my $self = shift;
    return (exists $self->{_hash}{SPID}) ?
	$self->{_hash}{SPID} :
	"na";
}
sub getNcbiOrgName {    
    my $self = shift;
    return (exists $self->{_hash}{NCBI_Organism_Name}) ?
	$self->{_hash}{NCBI_Organism_Name} :
	"na";
}
sub getPiEmail {    
    my $self = shift;
    return (exists $self->{_hash}{PI_EMAIL}) ?
	$self->{_hash}{PI_EMAIL} :
	"na";
}
sub getSpecies {    
    my $self = shift;
    return (exists $self->{_hash}{Species}) ?
	$self->{_hash}{Species} :
	"na";
}
sub getDraftSubmissionEnd {    
    my $self = shift;
    return (exists $self->{_hash}{Draft_Submission_End}) ?
	$self->{_hash}{Draft_Submission_End} :
	"na";
}
sub getEmbargoEndDate {    
    my $self = shift;
    return (exists $self->{_hash}{Embargo_End_Date}) ?
	$self->{_hash}{Embargo_End_Date} :
	"na";
}
sub getDraftSubmissionStart { 
    my $self = shift;
    return (exists $self->{_hash}{Draft_Submission_Start}) ?
	$self->{_hash}{Draft_Submission_Start} :
	"na";
}
sub getNcbiTaxId {    
    my $self = shift;
    return (exists $self->{_hash}{NCBI_Tax_ID}) ?
	$self->{_hash}{NCBI_Tax_ID} :
	"na";
}
sub getPiName {    
    my $self = shift;
    return (exists $self->{_hash}{PI_NAME}) ?
	$self->{_hash}{PI_NAME} :
	"na";
}
sub getProgramYear {    
    my $self = shift;
    return (exists $self->{_hash}{ProgramYear}) ?
	$self->{_hash}{ProgramYear} :
	"na";
}
sub getPmoProjId {    
    my $self = shift;
    return (exists $self->{_hash}{PMO_PROJECT_ID}) ?
	$self->{_hash}{PMO_PROJECT_ID} :
	"na";
}
sub getFinalSumbissionStart {    
	my $self = shift;
	return (exists $self->{_hash}{Final_Submission_Start} ) ?
	    $self->{_hash}{Final_Submission_Start} :
	    "na";
}
sub getGenus {
    my $self = shift;
    return (exists $self->{_hash}{Genus}) ?
	$self->{_hash}{Genus} :
	"na";
}
sub getPiHtml {    
    my $self = shift;
    return (exists $self->{_hash}{PI}) ?
	$self->{_hash}{PI} :
	"na";
}
sub getStrain {    
    my $self = shift;
    return (exists $self->{_hash}{Strain}) ?
	$self->{_hash}{Strain} :
	"na";
}
sub getOrganismName {    
    my $self = shift;
    return (exists $self->{_hash}{Name}) ?
	$self->{_hash}{Name} :
	"na";
}
sub getOrganismNameNoSpaces {    
    my $self = shift;
    my $oName = (exists $self->{_hash}{Name}) 
	?  $self->{_hash}{Name} 
        :  "na";
    $oName =~ s/ /_/g;
    return $oName;
}
sub getDataReleaseControl {    
    my $self = shift;
    return (exists $self->{_hash}{Data_Release_Control} ) ?
	$self->{_hash}{Data_Release_Control} :
	"na";
}
sub getAll {
   my $self = shift;
   return $self->{_hash};
}









1;


