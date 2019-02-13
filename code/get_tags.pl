#!/usr/bin/perl
use strict;
use Getopt::Long;

my $file    = '';
my $annot   = '';
my $indices = '';
my $cutoff  = 0;
my $help    = 0;

GetOptions ('files=s' => \$file,
	 					'annot=s' => \$annot,
						'indices=s' => \$indices,
						'cutoff=i' => \$cutoff,
						'help' => \$help);

########################################################
# USAGE
#
my $USAGE =<<USAGE;

     Usage:

         perl get_tags.pl [ -files -annot -indices -cutoff] [-help]

         where:

          files: directory that contains the fastq files (zipped or unzipped)
          annot: file that contains the lookup table of all possible tags
          indices: file that contains the sequences (8 bp) for each pGWDI index
          cutoff: cutoff for the number of reads; tags with less reads are discarded
          help: Prints out this helpful message

         If you are unsure how these files should look like go to REMIseq.org and download an example!

USAGE
#
######################################################
if($help){
    print "$USAGE\n";
    exit 0;
}

# DpnII: GATCCGTTGGA CTGCTG CGTGATTATGTATAATTT
# NlaIII:  CATGCGTTGGA CTGCTG CGTGATTATGTATAATTT

my($sample,$c,$orient,$single,$multiple,$inv,$not_found,$perc,$tag,$revcomp,$insert,$d,$entry,$st,$b,$all_reads,$DpnII,$NlaIII,$D,$count,$count_other,$line,$entry1,$entry2);
my(@div,@div1,@div2,@div3,@files,@help1);
my(%tags,%annots,%counts,%ind,%stats,%all);

#read in all .fastq files in current directory
opendir(DIR, $file) or die "cannot open directory";
my @files = readdir(DIR);
print "read files\n";

###############################################################
# read in position infos	                                    #
# possible upstream and reverse comlpement of downstream tags #
# are stored in hash %annots																  #
###############################################################

open(IN,$annot);
	while(<IN>){
		chomp;
		@div=split(/\t/,$_);
		$annots{uc $div[2]}{$div[0]."_".$div[1]}=$_;
		$revcomp = reverse uc $div[3];
		$revcomp =~ tr/ATGCatgc/TACGtacg/;
		$annots{uc $revcomp}{$div[0]."_".$div[1]}=$_;
	}
close(IN);

# count how often each taag occurs

foreach $entry (sort keys %annots){
	foreach $entry1 (sort keys %{$annots{$entry}}){
			if(! exists $counts{$entry}){
				$counts{$entry} = 1;
			} else {
				$counts{$entry}++;
			}
	}
}

# read indices for different pGWDI inserts
open(IN,$indices);
	while(<IN>){
		chomp;
		@div=split(/\t/,$_);
		$ind{$div[0]}=$div[2]."\t".$div[3];
	}
close(IN);
print "read annotations\n";

#######################################################
# go through each file and find reads containing      #
# end of vector (11 bp)	                              #
# tags and pGWDI index are extracted from those reads	#
#######################################################
$c = 0;
foreach $entry (@files){
	
	if (index($entry, ".fastq") != -1) {
		print $entry,"\n";
		undef(%tags);
		$all_reads=0;
		$count=0;
		$count_other=0;

		# read files either zipped or unzipped
		if ($entry =~ /.gz$/) {
	            open(IN, "gunzip -c $file/$entry |") || die "cannot open pipe to $entry";
	  	} else {
	            open(IN, "$file/$entry") || die "cannot open $entry";
	  	}

		while($line=<IN>){
				$line = <IN>;
				$all_reads++;
				chomp($line);
				
				$DpnII = index($line,"GATCCGTTGGA");
				$NlaIII  = index($line,"CATGCGTTGGA");

				# Does the read contain the end of the insert?
				if($DpnII > -1 | $NlaIII > -1){
					
					if($DpnII > -1){
						$st = $DpnII;
					}else{
						$st = $NlaIII;
					}
					$count++;
					# extract tag including GATC or CATG
					$b = substr($line,0,$st+4);
					$insert = substr($line,$st+17,6);
					# Does the tag have the correct size?
					if(length($b) >= 18 && length($b) <= 20 && exists $ind{$insert}){
						#print $line,"\n";
						# add insert sequence; add to hash; count how often it occurs
						$d = $b."_".$insert;
						if(!(exists $tags{$d})){
							$tags{$d} = 1;
						}else{
							$tags{$d} = $tags{$d} + 1;
						}
					}
				}elsif($line =~ /\A[acgt]+\z/i){
					$count_other++;
				}
				$line = <IN>;
				$line = <IN>;
			}
		close(IN);

		#number of reads
		@help1 = keys %tags;
		$D	   = @help1;
		$single = 0;
		$inv = 0;
		$multiple = 0;
		$not_found = 0;
		$sample = $entry;
		$sample =~ s/\.fastq\.gz//;
		$sample =~ s/\.fastq//;
		# add annotations
		foreach $tag (sort keys %tags){
			if($tags{$tag} > $cutoff){
				@div1 = split(/_/,$tag);
				@div3 = split(/\t/,$ind{$div1[1]});
				if($counts{$div1[0]} == 1){
					foreach $entry2 (sort keys %{$annots{$div1[0]}}){
						@div2 = split(/\t/,$annots{$div1[0]}{$entry2});
						$orient = get_orient(@div1,@div2,@div3);
						$all{$c} = $div1[0]."\t".$div1[1]."\t".$ind{$div1[1]}."\t".$tags{$tag}."\t".$annots{$div1[0]}{$entry2}."\tsingle\t".$orient."\t".$sample."\n";
						$single++;
					}
				} elsif($counts{$div1[0]} > 1){
					foreach $entry2 (sort keys %{$annots{$div1[0]}}){
						@div2 = split(/\t/,$annots{$div1[0]}{$entry2});
						$orient = get_orient(@div1,@div2,@div3);
						if($div2[0] eq "DDB0232429" && (($div2[1] > 2263132 && $div2[1] < 3015703) | ($div2[1] > 3016083 && $div2[1] < 3768654))){
							$all{$c} = $div1[0]."\t".$div1[1]."\t".$ind{$div1[1]}."\t".$tags{$tag}."\t".$annots{$div1[0]}{$entry2}."\tinverted repeat\t".$orient."\t".$sample."\n";
							$inv++;
						} else{
							$all{$c} = $div1[0]."\t".$div1[1]."\t".$ind{$div1[1]}."\t".$tags{$tag}."\t".$annots{$div1[0]}{$entry2}."\tmultiple\t".$orient."\t".$sample."\n";
							$multiple++
						}
					}
				}else{
					foreach $entry2 (sort keys %{$annots{$div1[0]}}){
						$all{$c} = $div1[0]."\t".$div1[1]."\t".$ind{$div1[1]}."\t".$tags{$tag}."\t\t\tno match"."\t".$sample."\n";
						$not_found++;
					}
				}
				$c++;
			}


		}

		# save stats for this run in hash %stats
		$perc = sprintf "%.2f", $count/$all_reads;
		$stats{$entry} = $all_reads."\t".$count."\t".$perc."\t".$single."\t" . $inv/2 ."\t".$multiple;


	}
}

# write to positions file (all)
open(OUT,">".$file."_positions_separate");
	foreach $entry (sort keys %all){
		print OUT $all{$entry};
	}
close(OUT);

# write stats file
open(OUT,">".$file."_stats");
	print OUT "filename","\t","reads","\t","reads_with_tags","\t","percentage_reads_with_tags","\t","unique_tags","\t","tags_in_inverted_repeat","\t","not_unique_tags","\n";
	foreach $entry (sort keys %stats){
		print OUT $entry,"\t",$stats{$entry},"\n"
	}
close(OUT);

########## subs ###############################################################

# @div1 : array containing [0] tag and [1] index
# @div2 : annotation for this tag
# @div3 : info about index (is it the right or left index?)
# returns orientation of insert

sub get_orient(@div1,@div2,@div3) {
	if($div3[1] eq "L"){
		if($div1[0] eq uc $div2[2]){
			$orient = "+";
		}else{
			$orient = "-";
		}
	}else{
		if($div1[0] eq uc $div2[2]){
			$orient = "-";
		}else{
			$orient = "+";
		}
	}
  return $orient;
}
