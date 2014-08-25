#!/usr/bin/perl

while (<>){
	my $old = <>;
	$old =~ s/ tag ".*?";//g;
	my $new = ($old =~ s/.*pfamA_acc "(.*?)".*gene_id "(.*?)".*exon_id "(.*?)"/$2_$3_$1/r);
	chomp($new);chomp($old);
	print $old," uniq_id ",$new,"\n";
}
