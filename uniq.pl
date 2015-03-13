#!/usr/bin/perl

while (<>){
	my $old = $_;
	$old =~ s/ tag ".*?";//g;
	my $new = ($old =~ s/.*pfamA_id "(.*?)".*pfamA_auto_reg "(.*?)".*gene_id "(.*?)".*exon_id "(.*?)"/$3_$4_$2_$1/r);
	my $new2 = "";
	if (/ccds_id ".*?";/) {$new2 = ($old =~ s/.*(ccds_id ".*?").*/$1\;/r);}
	$old =~ s/ ccds_id ".*?";*//g;
	$old =~ s/"ens_up_offset.*?";*//g;
	chomp($new);chomp($old);chomp($new2);
	if (!eof()) {if ($new2=="") {print $old," uniq_id ",$new," ",$new2,"\n";} else {print $old," uniq_id ",$new,"\n";}}
}