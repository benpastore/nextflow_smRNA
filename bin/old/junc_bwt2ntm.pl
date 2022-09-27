#!/usr/bin/perl

die "perl $0 <uniq.reads> <bowtie out>\n" if @ARGV<2;
my %nreads;
open ur,$ARGV[0];
while(<ur>)
{
	chomp;
	split;
	$nreads{$_[0]}=$_[1];
}
close ur;

my %nloc;
open in,$ARGV[1];
while(<in>)
{
	chomp;
	my ($seq)=(split)[0];
	$nloc{$seq}++;
}
close in;

open in,$ARGV[1];
while(<in>)
{
        chomp;
        my ($seq,$strand,$junc,$pos,$len)=(split)[0..4];
	my ($chr,$start,$junc_pos,$end)=(split/\|/,$junc)[0..3];
	my ($j_start,$j_ent)=split/-/,$junc_pos;
	print $chr,"\t",$start+$pos,"\t",$start+29,"\t",$seq,"\t",$nreads{$seq}/($nloc{$seq}*2),"\t",$strand,"\n";
	print $chr,"\t",$j_ent,"\t",($j_ent+$pos+length($len)-29),"\t",$seq,"\t",$nreads{$seq}/($nloc{$seq}*2),"\t",$strand,"\n";
}
close in;
