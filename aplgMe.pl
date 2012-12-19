#!/usr/bin/env perl
use strict;
use warnings;
use File::Temp;
use File::Copy;
use File::Path;
use File::Path qw(rmtree);
use File::Path qw(mkpath);
use File::Basename;
use Cwd;
use Cwd qw (realpath getcwd abs_path);
use lib "/house/homedirs/k/klabutti/lib";
use k;
use UGE;
use lib "/house/groupdirs/QAQC/scripts/versions/DEFAULT/lib";
use PGF::Utilities::Properties;
use PGF::Utilities::JobManager;
use Getopt::Long;
use vars qw/$optHelp $optDebug $optEstGenLen $optThreads
$optFrag $optFragInsert $optFragInsertStdDev $optAplgBin
$optJump $optJumpInsert $optJumpInsertStdDev $optName
$optInwardJump $optOutwardJump $optPloidy $optFragCov 
$optJumpCov $optEvalModeStandard $optEvalModeFull 
$optEvalModeCheat $optPhred64 $optPhred32  $optExtraOpts 
$optCleanUp $optReference $optWorkSpace $optIgnoreQual/;


#
# Kurt M. LaButti
#
#


#1.3
#modified qual format handling again so its better

#19nov2012  moved fail due to diff quals down so csvs would still be created
#16nov2012 v_1.2
#-added phred64 override, additional aplg args string
#-added cleanup

#15nov2012 modifications.   v_1.1
#-fixed $PATH export
#-added error handling if prep or run fail
#-updated default bin to /usr/common/jgi/assemblers/allpaths-lg/42816/bin
#-requirement to specify jump read orientation with -jin -jout
#-ability to provide a reference
#-eval mode switches.   
#    default=basic
#    if a reference is given default=standard
#    -full
#    -cheat



#19october2012 creation


my $cmdLine   = join " ", $0, @ARGV;

if( !GetOptions(
         "h"           => \$optHelp,
         "d"           => \$optDebug,
         "w=s"         => \$optWorkSpace,
         "b=s"         => \$optAplgBin,
         "l=i"         => \$optEstGenLen,
         "c"           => \$optCleanUp,


         "p32"         => \$optPhred32,
         "p64"         => \$optPhred32,
         "x=s"         => \$optExtraOpts,

         "r=s"        => \$optReference,
         "es"         => \$optEvalModeStandard,
         "ef"         => \$optEvalModeFull,
         "ec"         => \$optEvalModeCheat,

         "n=s"         => \$optName,
         "p=i"         => \$optPloidy,
         "t=i"         => \$optThreads,
         "f=s"         => \$optFrag,
         "fi=i"        => \$optFragInsert,
         "fs=i"        => \$optFragInsertStdDev,
         "fc=i"        => \$optFragCov,
         "jin"         => \$optInwardJump,
         "jout"        => \$optOutwardJump,
         "j=s"         => \$optJump,
         "ji=i"        => \$optJumpInsert,
         "js=i"        => \$optJumpInsertStdDev,
         "jc=i"        => \$optJumpCov,
         "iq"          => \$optIgnoreQual

    )
    ) { 
    usage("no options");
} 

if ($optHelp) { usage("help menu")};
unless ($optFrag) { usage("help menu")};
unless ($optJump) { usage("help menu")};
unless ($optInwardJump || $optOutwardJump) { usage("help menu")};
if ($optEvalModeStandard) {
    usage("you can only specify one evaulation mode!") 
	if ($optEvalModeFull);
    usage("you can only specify one evaulation mode!") 
	if ($optEvalModeCheat);
}
if ($optEvalModeFull) {
    usage("you can only specify one evaulation mode!") 
	if ($optEvalModeStandard);
    usage("you can only specify one evaulation mode!") 
	if ($optEvalModeCheat);
    
}
if ($optEvalModeCheat) {
    usage("you can only specify one evaulation mode!") 
	if ($optEvalModeFull);
    usage("you can only specify one evaulation mode!") 
	if ($optEvalModeStandard);
}
if ($optAplgBin) {
    usage("aplg bin: $optAplgBin location is bad!") unless (-d "$optAplgBin");
}
#==============================================================================#
# INITIALIZE INITIALIZE INITIALIZE INITIALIZE INITIALIZE INITIALIZE INITIALIZE #
#==============================================================================#
my $workSpace = ($optWorkSpace)
    ? realpath($optWorkSpace)
    : getcwd();
chomp $workSpace;
$optName =~ s/\s+/_/g if ($optName);

my $PRE = $workSpace;
my $REF = ($optName)
    ? $optName
    : 'aplgMe_assem';
my $DATA = 'allpaths';
my $RUN  = ($optName)
    ? "run_${optName}" 
    : 'run';
my $scriptsDir = "$PRE/$REF/$DATA/scripts";
my $frag     = (-s "$optFrag") 
    ? $optFrag
    : usage("bad frag: $!\n");
$frag = realpath($frag); chomp $frag;
my $fragIns  = ($optFragInsert)
    ? $optFragInsert
    : '250';
my $fragStdv = ($optFragInsertStdDev)
    ? $optFragInsertStdDev
    : '50';
my $fragCov = ($optFragCov)
    ? $optFragCov
    : 50;
my $jump     = (-s "$optJump") 
    ? $optJump
    : usage("bad frag: $!\n");
$jump = realpath($jump); chomp $jump;
my $jumpIns  = ($optJumpInsert)
    ? $optJumpInsert
    : 4000;
my $jumpStdv = ($optJumpInsertStdDev)
    ? $optJumpInsertStdDev
    : 500;
my $jumpCov = ($optJumpCov)
    ? $optJumpCov
    : 50;
my $jumpDir = ($optInwardJump)
    ? 'inward'
    : ($optOutwardJump)
    ?'outward'
    : "outward";
my $aplgBin = ($optAplgBin) 
    ? $optAplgBin
    : "/usr/common/jgi/assemblers/allpaths-lg/42816/bin";
my $prep_script = (-s "$aplgBin/PrepareAllPathsInputs.pl")
    ?  "$aplgBin/PrepareAllPathsInputs.pl"
    : usage("prep script failure");
my $run_script = (-s "$aplgBin/RunAllPathsLG")
    ?  "$aplgBin/RunAllPathsLG"
    : usage("run script failure");
my $fasta2fastb_script = (-s "$aplgBin/Fasta2Fastb")
    ? "$aplgBin/Fasta2Fastb"
    : usage("fasta2fasb script failure");
my $estGenomeLen = ($optEstGenLen) 
    ? $optEstGenLen
    : 50000;
my $ploidy = ($optPloidy)
    ? $optPloidy
    : 1;
my $threads = ($optThreads)
    ? $optThreads
    : 8;
my $picard    = "/jgi/tools/misc_bio/picard/versions/picard-tools-1.48";
my $reference = ($optReference)
    ? realpath($optReference)
    : '';
chomp $reference;
my $evalMode  = ($optEvalModeStandard) ? "STANDARD"
    : ($optEvalModeFull)  ? "FULL"
    : ($optEvalModeCheat) ? "CHEAT"
    : "BASIC";  #default, none specified
#if no eval was specified, but there is a ref, do standard by default
$evalMode = ($evalMode =~ /BASIC/ && $reference)
    ? "STANDARD"
    : "BASIC";
#==============================================================================#
# VALIDATE VALIDATE VALIDATE VALIDATE VALIDATE VALIDATE VALIDATE VALIDATE VALI #
#==============================================================================#
#exit if trying to use standard, full, or cheat w/o a reference
if ($optEvalModeStandard || $optEvalModeFull || $optEvalModeCheat) {
    usage("-r=s required if using evaluation mode other than BASIC!") 
	unless ($optReference);
}
#check ref is supplied
if ($optReference) {
    usage("reference [$reference] does not exist!") 
	unless (-s "$reference");
}
#==============================================================================#
#  MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN  #
#==============================================================================#
k::msg("running: $cmdLine");
k::msg("evalMode=$evalMode");
k::msg("reference=$reference") if ($optReference);
k::msg("workSpace=$workSpace");
unless (-d $workSpace) { 
    mkdir($workSpace,0777) or die "Can't make workdir: $workSpace: $!\n";
}

#setup dir
my $setup = setupAPLG($scriptsDir,$REF,$workSpace,$DATA);
if ($setup == 0) {k::msg("setup failure"); exit;}

chdir "$PRE/$REF/$DATA";
my $cwd = getcwd();
k::msg("CURR_DIR=$cwd");

#process reference if supplied 
processRef($reference) if ($optReference);

#create csv
my ($groups, $libs, $phred64) = createCsv($frag,$jump);

#prep and run
runAssem($groups, $libs, $phred64);

#cleanup
if ($optCleanUp) {
rmtree "$cwd/read_cache";
unlink glob "frag_reads*";
unlink glob "jump_reads*";
}
#==============================================================================#
# soubroutines soubroutines soubroutines soubroutines soubroutines soubroutines#
#==============================================================================#
sub setupAPLG {
    my ($scriptsDir, $REF, $assemSpace, $DATA) = @_;
    my $script = "/global/u1/k/klabutti/bin/setupAPLG.pl";
    my $cmd = "$script "
	."-s $scriptsDir "
	."-bio \"$REF\" "
	."-space $assemSpace "
	."-data $DATA ";
    
    k::msg("setupAPLG: $cmd");
    system($cmd);
    
    #
    # check if dir created
    unless (-d "$assemSpace/$REF/$DATA") {
	k::msg("setupAPLG FAILED!: [$assemSpace/$REF/$DATA]");
	return "0";
    } else {
	return "1";
    }
}
#==============================================================================#
sub processRef {
    my $ref = shift @_;
    
#copy to REF dir
    k::msg("copying reference to: $PRE/$REF/genome.fasta");
    copy("$ref", "$PRE/$REF/genome.fasta");

#convert to fastb in REF dir
    my $cmd = "$fasta2fastb_script IN=$ref OUT=$PRE/$REF/genome.fastb";
    k::msg("processing refernce: $cmd");
    system($cmd);
}

#==============================================================================#
sub createCsv {
    my ($frag,$jump) = @_;
    my $phred64;
    
    #check qual format
    my $fq = k::guessFastqFormat($frag);
    my $jq = k::guessFastqFormat($jump);
  

    my $libs = "library_name, "
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

    $libs .= 
	"frag,aplg,${REF},fragment,1,${fragIns},${fragStdv},,,inward,0,0\n";
    $libs .= 	
	"jump,aplg,${REF},jumping,1,,,${jumpIns},${jumpStdv},$jumpDir,0,0\n";

    
    my $groups = "group_name, library_name, file_name\n";
    $groups .= "frag,frag,$frag\n";
    $groups .= "jump,jump,$jump\n";


    open G, ">$PRE/$REF/$DATA/in_groups.csv" or die "Can't open: $!\n";
    print G "$groups";
    open L, ">$PRE/$REF/$DATA/in_libs.csv" or die "Can't open: $!\n";
    print L "$libs";
    close G;
    close L;

    exit "in_groups.csv issue" 
	unless (-s "$PRE/$REF/$DATA/in_groups.csv");
    exit "in_libs.csv issue" 
	unless (-s "$PRE/$REF/$DATA/in_libs.csv");
    


#auto handle quals
    if ($fq =~ /sanger/ && $jq =~ /sanger/) {
	$phred64 = 'False';
    } elsif ($fq =~ /ill/ && $jq =~ /ill/) {
	$phred64 = 'True';
    } else {
	
	if ($optPhred64 || $optPhred32) {
	    #override switches
	    $phred64 = "True" if $optPhred64;
	    $phred64 = "False" if $optPhred32;
	} else { 
	    usage("MIX OF QUALITY FORMATS FOUND frag=[$fq] jump=[$jq];  try using setting the phred offset manually with -jp32, -fp32, -jp64 or , -fp64");
	}
    }
    
    return "$PRE/$REF/$DATA/in_groups.csv", 
    "$PRE/$REF/$DATA/in_libs.csv",
    $phred64;
}

#==============================================================================#
sub runAssem {
    my ($groups, $libs, $phred64) = @_;

    my $pathCmd = "export PATH=$aplgBin/:\$PATH"; 
    
    my $prepCmd = "$prep_script "
	."PICARD_TOOLS_BIN=$picard "
	."DATA_DIR=$PRE/$REF/$DATA "
	."IN_GROUPS_CSV=$groups "
	."IN_LIBS_CSV=$libs "
	."GENOME_SIZE=$estGenomeLen "
	."FRAG_COVERAGE=$fragCov "
	."JUMP_COVERAGE=$jumpCov "
	."PHRED_64=$phred64 "
	."PLOIDY=$ploidy ";
    k::msg("prep=[$prepCmd]");
    open  O, ">$cwd/prep.sh" or die "Can't open prep.sh: $!\n";
    print O "$pathCmd; $prepCmd > prep.log";
    close O;
    chmod 0777, "$cwd/prep.sh";
    system("$pathCmd; $prepCmd > prep.log");    

    checkSuccess("$PRE/$REF/$DATA/frag_reads_orig.fastb");
    checkSuccess("$PRE/$REF/$DATA/jump_reads_orig.fastb");

    my $runCmd = "$run_script "
	."PRE=$PRE "
	."REFERENCE_NAME=$REF "
	."DATA_SUBDIR=$DATA "
	."RUN=$RUN "
	."THREADS=$threads "
	."VAPI_WARN_ONLY=True "
	."OVERWRITE=True";

	$runCmd .= " EVALUATION=$evalMode" if ($optReference);
	$runCmd .= " $optExtraOpts" if ($optExtraOpts);

    k::msg("run=[$runCmd]");
    open  OR, ">$cwd/run.sh" or die "Can't open run.sh: $!\n";
    print OR "$pathCmd; $runCmd > run.log";
    close OR;
    chmod 0777, "$cwd/run.sh";
    system("$pathCmd; $runCmd > run.log");

    checkSuccess("$PRE/$REF/$DATA/$RUN/ASSEMBLIES/test/assembly.report");

}

#==============================================================================#
sub checkSuccess {
    my $expectedFile = shift @_;

    unless (-s $expectedFile) {
	warn "expected file does not exist; [$expectedFile], exiting!\n";
	exit;
    }
}

#==============================================================================#
sub usage {
    my ($msg) = @_;
 
    unless ($msg) {$msg = "";}
 
print <<USE;
 
$msg

*version 1.3*
 
$0 [options] -f <fragFastq> -j <jumpFastq> 
   -h  help menu
   -d  debug

   -f=s    frag fastq         REQUIRED
   -j=s    jump fastq         REQUIRED
   -jin OR -jout    jump ori  REQUIRED
 
   #basic options
   -b=s  path to aplg bin  (def=/usr/common/jgi/assemblers/allpaths-lg/42816/bin)
   -l=i  estimated genome len (def=50000)
   -n    assembly proj name   (def=APLG_ASSEM)
   -t    num threads for aplg (def=8)
   -x=s  extra options string for AllPaths.  Enclose in " "   
   -p64  auto qual guess override (phred_64=True) 
   -p32  auto qual guess overried (phred_64=False)
   -c    clean up  read_cache and prepped read data  

   #fragment options
   -fi=i   frag insert size   (def=250)
   -fs=i   frag insert stdDev (def=50)
   -fc=i   frag coverage      (def=50)

   #jumping options
   -ji=i   jump insert size   (def=4000)
   -js=i   jump insert stdDev (def=500)
   -jc=i   jump coverage      (def=50)

   #reference
   -r=s    reference fasta (def=no reference)

   #evaultion modes, REQUIRES -r, default=basic with no reference
      *** use ONE of the following if specifying a mode ***
    -es    STANDARD evaluation; run eval modules using supplied reference
    -ef     FULL evaluation; turn on in-place evaluation in certain modules
    -ec     CHEAT evaluation; use reference to guide the assembly slightly


 Purpose:  a lightweight wrapper for running the APLG assembler.
           By default geared toward assembling mito, but can be 
           changed with options.  Runs command line (no cluster).

 
USE
exit;
} 

