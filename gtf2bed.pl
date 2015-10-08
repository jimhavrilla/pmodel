#!/usr/bin/perl

while (<>){
	my $old = $_;
	$chr=($_=~s/^(\w*)\s*(\d*)\s*(\d*)\s*(.*\n)/$1/r);
	$start=($_=~s/^(\w*)\s*(\d*)\s*(\d*)\s*(.*\n)/$2/r);
	$end=($_=~s/^(\w*)\s*(\d*)\s*(\d*)\s*(.*\n)/$3/r);
	$new=($_=~s/^(\w*)\s*(\d*)\s*(\d*)\s*(.*\n)/$4/r);
	if (!eof()) {print $chr,"\t",($start-1),"\t",($end),"\t",$new;}
}
