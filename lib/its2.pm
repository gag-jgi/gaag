package its2;

=head1 NAME

PMO

=head1 SYNOPSIS

Ex usage:
#!/jgi/tools/bin/perl -w
use strict;
use lib "/home/klabutti/bin/modules";
use ITS;
 
my $projId = shift @ARGV;

my $itsObj = ITS->new("$projId");

my $piName = $itsObj->getPiName();

print "PI name: $piName\n";


The available functions (pretty self explanatory):

sub getPiEmail
sub getPiName
sub getPiId 
sub getPropId
sub getPropTitle
sub getAll


=head1 DESCRIPTION
This module accesses ITS database for collab. info given a project id

=head1 AUTHOR(s)

Kurt M. LaButti

=head1 HISTORY

26sept2012; creation

=over

=item *

Kurt M. LaButti 26September2012 Creation

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

    $self->_getIts("$projectId");
    
    return $self;

}

#============================================================================#
sub _getIts {
    
    my $self = shift;    
    my $seqProjId = shift;
    my %info;
    
    $ENV{ORACLE_HOME} = '/jgi/tools/oracle_client/DEFAULT';
    $ENV{TNS_ADMIN} = '/jgi/tools/oracle_client/DEFAULT/network/admin';    
    
    my $dbh;
    eval {
	$dbh = DBI->connect( 'dbi:Oracle:dwprd1', 'dw_user', 'dw_user' );
    };
    if ($@) {
	die $@. "\n";
    }
    
    
    my $sql = "select pr.proposal_id, pr.title, pr.principal_investigator_name, pr.principal_investigator_id, sq.taxonomy_info_id, sq.sequencing_project_id, sq.sequencing_project_name
from 
dw.proposal pr, 
dw.sequencing_project sq,
dw.final_deliv_project fd
where 
sq.sequencing_project_id =  $seqProjId and
pr.proposal_id = fd.proposal_id and 
fd.proposal_id = pr.proposal_id and
sq.final_deliv_project_id = fd.final_deliv_project_id";
    
    
    my $sth = $dbh->prepare("$sql");
    $sth->execute()
	or confess "Couldn't execute statement: " . $sth->errstr;
#    my $numRows = $sth->rows();
#    print "numRows: $numRows\n";
    

    my ($propId,$title,$inv,$invId,$taxId,$proj_id,$proj_name) = $sth->fetchrow_array;
    

    confess "No DATA FOUND!\n" unless ($title);
    
#get email from contact table
    $sql = "select email_address from contact where contact_id = $invId";
    $sth = $dbh->prepare("$sql");
    $sth->execute()
	or confess "Couldn't execute statement: " . $sth->errstr;
    my ($email) = $sth->fetchrow_array;
    
#get ncbi_name from taxonomy info table
 
$sql = "select ncbi_organism_name from taxonomy_info where taxonomy_info_id = $taxId";
$sth = $dbh->prepare("$sql");
$sth->execute()
    or die "Couldn't execute statement: " . $sth->errstr;
my ($name) = $sth->fetchrow_array;
print "name: $name\n";


#get final_deliv_name
$sql = "select ncbi_organism_name from taxonomy_info where taxonomy_info_id = $taxId";



    
    $sth->finish();
    $dbh->disconnect();
    
#stuff into hash
    $info{pi_name} = $inv;
    $info{pi_email} = $email;
    $info{pi_id} = $invId;
    $info{proposal_title} = $title;
    $info{proposal_id} = $propId;
    $info{taxonomy_id} = $taxId;
    $info{ncbi_name}   = $name;
    $info{proj_name}   = $proj_name;
    $self->{_hash} = \%info;
}



#============================================================================#

sub getPiEmail {    
    my $self = shift;
    return (exists $self->{_hash}{pi_email}) ?
	$self->{_hash}{pi_email} :
	"na";
}

sub getPiName {    
    my $self = shift;
    return (exists $self->{_hash}{pi_name}) ?
	$self->{_hash}{pi_name} :
	"na";
}

sub getPiId {    
    my $self = shift;
    return (exists $self->{_hash}{pi_id}) ?
	$self->{_hash}{pi_id} :
	"na";
}

sub getPropId {    
    my $self = shift;
    return (exists $self->{_hash}{proposal_id}) ?
	$self->{_hash}{proposal_id} :
	"na";
}

sub getPropTitle {    
    my $self = shift;
    return (exists $self->{_hash}{proposal_title}) ?
	$self->{_hash}{proposal_title} :
	"na";
}

sub getTaxId {    
    my $self = shift;
    return (exists $self->{_hash}{taxonomy_id}) ?
	$self->{_hash}{taxonomy_id} :
	"na";
}

sub getBioName {    
    my $self = shift;
    return (exists $self->{_hash}{ncbi_name}) ?
	$self->{_hash}{ncbi_name} :
	"na";
}

sub getAll {
    my $self = shift;
    return $self->{_hash};
}




1;
