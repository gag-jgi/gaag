package calc_cvg
use strict;
use warnings;

##
#KML 27feb08
##
#function for calculating number of reads needed for desired depth given 
#some basic information about a genome and its data set
#  (cvg * gen_size) / (Q20read_length * pass_rate) = num_reads 


sub calc_num_reads {
    my ($genome_size, $desired_cvg, $avg_read_len, $pass_rate) = @_;
    my $num_reads = ($desired_cvg * $genome_size) / ($avg_read_len * $pass_rate);
    return $num_reads;
}


sub calc_pass_rate {
 my ($genome_size, $cvg, $avg_read_len, $num_reads) = @_;
 my $pass_rate = ($cvg * $genome_size) / ($num_reads * $avg_read_len);
 return $pass_rate;
}


sub calc_cvg {
    my ($genome_size, $avg_read_len, $pass_rate, $num_reads) = @_;
    my $cvg = ($avg_read_len * $pass_rate * $num_reads) / $genome_size;
    return $cvg;
}




1;
