#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use lib '/Users/remo/src/BioPerl-1.6.1';
use Bio::SeqIO;
use Bio::SearchIO;

my $blastbinfold = '/Users/remo/src/ncbi-blast-2.2.25+/bin';
my $fasta = 'cab.fa';

open(OUT,">$fasta\_longest_stem.txt");
print OUT join("\t",qw(seqid longest qstart qend sstart ssend));
print OUT "\n";

my $seqio = Bio::SeqIO->new(-file => $fasta,
                             -format => 'fasta');

my $c = 1;

while(my $seq = $seqio->next_seq) {
  analyze_seq($seq,$c);
  $c++;
}

sub analyze_seq {
  my $seq = shift;
  my $c = shift;
  my $seqid = $seq->id;
  my $seqcount = "$c\.fa";
  my $blastout =  "$c\.blastn";
  my $log = "$c\.log";

  my $forio = Bio::SeqIO->new(-file => ">$seqcount",
                              -format => 'fasta');

  $forio->write_seq($seq);

  system("$blastbinfold\/makeblastdb -dbtype nucl -in $seqcount > $log");
  die "MAKEBLASTDB ERROR ON SEQUENCE $seqid\: $!" unless $? == 0;
  #system("$blastbinfold\/blastn -query $seqcount -db $seqcount -task blastn -dust no -word_size 6 -out $blastout -gapopen 2 -gapextend 1 -penalty -1 -reward 1 > $log");
  #die "BLASTN ERROR ON SEQUENCE $seqid\: $!" unless $? == 0;
  system("$blastbinfold\/blastn -query $seqcount -db $seqcount -task blastn -dust no -word_size 6 -out $blastout > $log");
  die "BLASTN ERROR ON SEQUENCE $seqid\: $!" unless $? == 0;

  parse_blastout($blastout,$seq->id);

  system("rm $seqcount\*");
  die "RM ERROR WHILE WORKING ON SEQUENCE $seqid\: $!" unless $? == 0;
  system("rm $blastout");
  die "RM ERROR WHILE WORKING ON SEQUENCE $seqid\: $!" unless $? == 0;
  system("rm $log");
  die "RM ERROR WHILE WORKING ON SEQUENCE $seqid\: $!" unless $? == 0;
}

sub parse_blastout {
  my $out = shift;
  my $seqid = shift;

  my $longestlength = 0;
  my $longest = 'NA';

  my $in = Bio::SearchIO->new(-format => 'blast', 
                              -file => $out);

  while(my $result = $in->next_result) {
    while(my $hit = $result->next_hit) {
      while(my $hsp = $hit->next_hsp) {

        next if $hsp->subject->strand != -1;

        my $qs = $hsp->query->start;
        my $qe = $hsp->query->end;
        my $ss = $hsp->subject->start;
        my $se = $hsp->subject->end;

        next if ($qs <= $se && $qe >= $ss);
        next if ($ss <= $qe && $se >= $qs);

        if($hsp->length > $longestlength) {
          $longest = $hsp;
          $longestlength = $hsp->length;
        }
      }
    }
  }

  if($longest eq 'NA') {
    print OUT join("\t",$seqid,0,0,0,0,0);
  }

  else {
    print OUT join("\t",$longest->query->seq_id,$longest->length('query'),$longest->query->start,$longest->query->end,$longest->subject->start,$longest->subject->end);
    print_hsp($longest);
  }

  print OUT "\n";
}

sub print_hsp {
  my $hsp = shift;
  print $hsp->query->seq_id."\t";
  print "HSP QUERY STRAND: ".$hsp->strand('query')."\t";
  print "HSP SUBJECT STRAND: ".$hsp->strand('subject')."\t";
  print "HSP LENGTH: ".$hsp->length."\t";
  print "QUERY:".$hsp->query->start."-".$hsp->query->end."\t";
  print "SUBJECT:".$hsp->subject->start."-".$hsp->subject->end."\n";
}
