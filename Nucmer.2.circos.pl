#!/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(floor ceil);
use Tie::IxHash;


#4/12/2015 - added code to order contigs by total length of links, rather than number of links - 
# (--ref_order='relw'; --query_order='relw' options).

my $min_hit_len = 1000;
my $min_ID = 95;
my $IS_coords = 0;
my $min_chr_len = 1000;
my $parse_contig_names = 0;
my $ribbons = 0;
my $promer = 0;
my $colour_links_by_query = 0;
my $debug = 0;
my $ref_order = 0;
my $query_order = 0;
my $no_query_labels = 0;
my $no_ref_labels = 0;
my $help = 0;
my $gradual_shading = 0;
my $col_start_h = 0;
my $col_start_s = 0.2;
my $col_end_s = 1.0;
my $col_start_v = 0.6;
my $col_end_v = 0.9;
my $default_v_fornongrad =0.95;
my $default_s_fornongrad = 0.95;
my $zoom_query = 0;
my $prefix = "";
my $label_size_param = "24p";

my $min_ref_chr_len = -1;
my $min_query_chr_len = -1;

my $colour_rainbow_links =-1;

my $flipquery= 1;
my $flipreference = 0;

my $ref_contigs_to_keep = "";
my $query_contigs_to_keep = "";

my $opt = GetOptions("prefix=s"=>\$prefix,"min_hit_len=i"=>\$min_hit_len,"min_ID=i"=>\$min_ID,"coord-file"=>\$IS_coords,
"min_chr_len=i"=>\$min_chr_len,"parse_contig_names"=>\$parse_contig_names,"ribbons"=>\$ribbons,"promer"=>\$promer,
"colour_links_by_query"=>\$colour_links_by_query,"debug"=>\$debug,"ref_order=s"=>\$ref_order,"query_order=s"=>\$query_order,
"no_query_labels"=>\$no_query_labels,"no_ref_labels"=>\$no_ref_labels,"help|?"=>\$help,"gradual_shading=i"=>\$gradual_shading,
"col_start_s=f"=>\$col_start_s,"col_end_s=f"=>\$col_end_s,"col_start_v=f"=>\$col_start_v,"col_end_v=f"=>\$col_end_v,
"col_start_h=f"=>\$col_start_h,"zoom_query=f"=>\$zoom_query,"label_size=s"=>\$label_size_param,"min_ref_chr_len=i"=>\$min_ref_chr_len,"min_query_chr_len=i"=>\$min_query_chr_len,
"colour_rainbow_links=n"=>\$colour_rainbow_links,"flipreference!"=>\$flipreference,"flipquery!"=>\$flipquery,"query_contigs_to_keep=s"=>\$query_contigs_to_keep,"ref_contigs_to_keep=s"=>\$ref_contigs_to_keep);



if ($help) { $opt =0 ; } 
#if ($ref_order || $query_order) { die("SORRY!!! - ref_order and query_order other than alphabetical don't work yet\n"); }
if (lc($ref_order) eq "rel" && lc($query_order) eq "rel") { print STDERR "ERROR! - Cannot specify 'rel' for both ref and query ordering - at least one needs to be determined\n\n"; $opt = 0; }

if (!$opt || scalar @ARGV != 1) { 
	print STDERR "usage is: perl $0 [options] <deltafile>\n\n";
	print STDERR "This script might/should/sometimes writes all the input files circos needs, but doesn't actually RUN circos. After running this, you will have a circos.conf and karyotype.txt and links.txt (maybe more) files in this folder. Then run circos (e.g. perl /nfs/users/nfs_j/jc17/software/circos-0.66/bin/circos) to generate cricos.png and circos.svg output files. This isn't just lazy - you might want to change circos parameters!\n\n";
	print STDERR "options are:\n";
	print STDERR "\n";
	print STDERR "INPUT OPTIONS:\n";
	print STDERR "\t--coord-file\tpass a coordinate file. I assume its -dtlro format - i.e. 13 columns\n";
	print STDERR "\t--promer\tdeltafile is from promer - expect extra columns!\n";
    print STDERR "\t--parse_contig_names\tby default, I try and keep the contig ID passed to circos the same as the ID in the input files.. but if you have special characters (e.g. ':') this upsets circos. Use this flag! - your printed names will be consistent, but your circos files won't use these, so it will  be harder to make sense of them / modify thing downstream of this script. Might make this default at some point in future..\n";	
	print STDERR "\t--query_contigs_to_keep\tby default, I keep all contigs that have links. If this is a filename, read this and just keep these\n";
	print STDERR "\t--ref_contigs_to_keep\tby default, I keep all contigs that have links. If this is a filename, read this and just keep these\n";
	
	print STDERR "\n";
	print STDERR "DISPLAY/OUTPUT OPTIONS:\n";
	print STDERR "\t--min_hit_len=1000\tlength of smallest hit (passed to show-coords with -L flag)\n";
	print STDERR "\t--min_chr_len=1000\tlength of smallest contig/scaffold/chromosome to include\n";
	print STDERR "\t--min_query_chr_len=1000\tlength of smallest QUERY contig/scaffold/chromosome to include (overrides min_chr_len)\n";
	print STDERR "\t--min_ref_chr_len=1000\tlength of smallest REF contig/scaffold/chromosome to include (overrides min_chr_len)\n";
	print STDERR "\t--ribbons\tuse ribbons rather than thin links\n";
	print STDERR "\t--gradual_shading=i\texperimental feature! - use i different colour intensities along scaffolds (default = 1)\n";
	print STDERR "\t--col_start_s=f\tstart 'saturation' level for gradual shading - also used as fixed value without gradual_shading (default ".$col_start_s." with gradual_shading, ".$default_s_fornongrad." without)\n";
	print STDERR "\t--col_end_s=f\tend 'saturation' level for gradual shading - default ".$col_end_s.", ignored otherwise\n";
	print STDERR "\t--col_start_v=f\tstart 'value' level for gradual shading, or fixed value without gradual shading - default ".$col_start_v."\n";
	print STDERR "\t--col_end_v=f\tend 'value' level for gradual shading - default ".$col_end_v."\n";
	print STDERR "\t--col_start_h=f\tstart 'hue' level for colour picker\n";
	print STDERR "\t--colour_links_by_query\tlinks have a different colour for each query scaffold.. default is for reference scaffold\n";
	print STDERR "\t--colour_rainbow_links=i\tcolour each link differently, as possible with i colours\n";
	print STDERR "\t--ref_order=file\teither a text file with identifiers in the order you want them to appear, or a fasta file with identifiers in the order you want them, otherwise sorted alphabetically, or relative to the query order (to try and minimise crossing - 'rel'; 'rel' is number of links, 'relw' is total length of links - probably better)\n";
	print STDERR "\t--query_order=file\tas above, but for query scaffolds\n";
	print STDERR "\t--flipquery\tflip query scaffolds so ribbons are untwisted etc. (default=yes, unless flipreference is used. To disable, use --noflipquery)\n";
	print STDERR "\t--flipreference\tflip reference scaffolds so ribbons are untwisted etc. (default=no, to disable, use --noflipquery)\n";
	print STDERR "\t--prefix=string\tadd this string to all output files\n";
	print STDERR "\t--zoom_query=f\tzoom query karyotype by this factor\n";
	print STDERR "\n";
	print STDERR "LABELLING OPTIONS:\n";
	print STDERR "\t--no_query_labels\tdon't show labels for query chromosomes\n";
	print STDERR "\t--no_ref_labels\tdon't show  labels for ref chromosomes\n";
	print STDERR "\t--label_size=s\tset size of labels. Include UNITS - i.e. 24p, not just 24 (default)\n";
	print STDERR "\n";
	print STDERR "MISCELLANEOUS OPTIONS:\n";
	print STDERR "\t--help|-h|-?\tprint this\b";
	print STDERR "\t--debug\tprint some debugging information\n";
	exit;
}

if ($flipreference) { $flipquery = 0; } 
if ($min_query_chr_len < 0 ) { $min_query_chr_len = $min_chr_len; }
if ($min_ref_chr_len < 0 ) { $min_ref_chr_len = $min_chr_len; }
if ($gradual_shading <= 1) { $col_start_s = $default_s_fornongrad; $col_start_v = $default_v_fornongrad;}


#########
#
# GENERATE A MUMMER 'COORDS' FILE IF NOT PASSED ON COMMANDLINE
#
##
my $coordsfile = "";
if (!$IS_coords) { 
$coordsfile = $ARGV[0].".coords";
my $commandline = "show-coords -dTlro -I ".$min_ID." -L ".$min_hit_len." ".$ARGV[0]." > ".$ARGV[0].".coords\n";
print `$commandline`;
}
else { $coordsfile = $ARGV[0]; } 

####
#
# PARSE COORDS FILE. BY DEFAULT (CURRENTLY ALWAYS, I ONLY PRINT  SCAFFOLDS THAT HAVE SOME LINKS.. THIS WILL PROBABLY NEED CHANGING AT SOME POINT
#
# MIGHT (parse_contig_names) NEED TO NOT USE ACTUAL INPUT FILE NAMES AS IDENTIFIERS (can't have some special chars), BUT GENERATE THEM AND
# STORE BACK-AND-FORTH MAPPING
#
# PRINT CIRCOS links.txt file as I go along
#
######

#g1 is reference
#g2 is query as passed to mummer
#hashes for contig lengths of genomes 1 and 2
my %contigs2len_g1;
my %contigs2len_g2;
my $q_scaff_count = 1;
my $r_scaff_count = 1;
#original name -> translation
my %q_scaff_translation;			
my %r_scaff_translation;
#translation -> original name
my %q_scaff_translation_back;
my %r_scaff_translation_back;


#contig 'keep list' for reference
my %r_scaff_keep;
#contig 'keep list' for query
my %q_scaff_keep;

if ($ref_contigs_to_keep ne "" ) { 
	open(RCIN,"<$ref_contigs_to_keep") or die("cannot read file ".$ref_contigs_to_keep);
	while (<RCIN>) { 
		chomp;
		$r_scaff_keep{$_}++;
	
	}
	close RCIN;
	print STDERR "read a total of ".(scalar keys %r_scaff_keep)." reference contigs to include from ".$ref_contigs_to_keep;
	if (scalar keys %r_scaff_keep < 1) { die("keeping no reference contigs.. doesn't make sense to continue under these circumstances"); }

}

if ($query_contigs_to_keep ne "" ) { 
	open(QCIN,"<$query_contigs_to_keep") or die("cannot read file ".$query_contigs_to_keep);
	while (<QCIN>) { 
		chomp;
		$q_scaff_keep{$_}++;
	
	}
	close QCIN;
	print STDERR "read a total of ".(scalar keys %q_scaff_keep)." query contigs to include from ".$query_contigs_to_keep;
	if (scalar keys %q_scaff_keep < 1) { die("keeping no query contigs.. doesn't make sense to continue under these circumstances"); }
}

## stores all links count / weights
# hash of hashes, ref->query->count
my %links;
# hash of hashes, query->ref->count
my %rlinks;

## stores all links weights!
my %links_weight;
my %rlinks_weight;

#swap weight
#are links forward or backward on each query scaffold!
#one number per query scaffold
#
my %links_swap_weight;

my $rscaff="NA";
my $qscaff="NA";
open(LI,">",$prefix."links.txt");
open(CO,"<$coordsfile") or die("cannot find coords file - does deltafile ".$ARGV[0]." exist? is nucmer properly installed\n");

my $correct_len = 12;
my $rname = 11;
my $qname = 12;
my $rlen = 7;
my $qlen = 8;

if ($promer) { 
	$correct_len = 14;
	$rname = 13;
	$qname = 14;
	$rlen = 9;
	$qlen = 10;
}
my $link_count = 0;
while (<CO>) { 
	chomp;
	my @line = split /\t/;

	#check i'm in a data line
	if  ( (scalar @line >= $correct_len ) && ( $line[0] =~ m/^[0-9]+/ ) ) { 
		if ($parse_contig_names) { 
			if (exists $r_scaff_translation{$line[$rname]} ) { 
				$rscaff = $r_scaff_translation{$line[$rname]};
			} else { 
				$rscaff = "R".$r_scaff_count;
				$r_scaff_translation{$line[$rname]} = $rscaff;
				$r_scaff_translation_back{$rscaff} = $line[$rname];
				$r_scaff_count++;
			}
			if (exists $q_scaff_translation{$line[$qname]} ) { 
				$qscaff = $q_scaff_translation{$line[$qname]};
			} else { 
				$qscaff = "Q".$q_scaff_count;		
				$q_scaff_translation{$line[$qname]} = $qscaff;
				$q_scaff_translation_back{$qscaff} = $line[$qname];
				$q_scaff_count++;
			}
		} else  {
			$rscaff = $line[$rname];
			$qscaff = $line[$qname];
		}
		
		my $ref_keep = 1;
		my $query_keep  = 1;
		
		if ( $query_contigs_to_keep ne "" ) { 
			if ( ! exists $q_scaff_keep{$line[$qname]} ) { $query_keep = 0; }
		}
		if ( $ref_contigs_to_keep ne "" )  {
			if ( ! exists $r_scaff_keep{$line[$rname]} ) { $ref_keep = 0; }
		}
		if ($line[$rlen] > $min_ref_chr_len && $ref_keep) { 
		#	$line[11] =~ s/\:/\|/;
			$contigs2len_g1{$rscaff} = $line[$rlen];
		}
		if ($line[$qlen] > $min_query_chr_len && $query_keep) {
		#	$line[12] =~ s/\:/\|/; 
			$contigs2len_g2{$qscaff} = $line[$qlen];
		}
		if ($line[$rlen] > $min_ref_chr_len && $line[$qlen] > $min_query_chr_len && $ref_keep && $query_keep) { 
			if ($colour_rainbow_links > 0) { 
				my $link_val = $link_count % $colour_rainbow_links;
				print LI $rscaff."\t".$line[0]."\t".$line[1]."\t".$qscaff."\t".$line[2]."\t".$line[3]." value=".$link_val."\n";
			} else { 
				print LI $rscaff."\t".$line[0]."\t".$line[1]."\t".$qscaff."\t".$line[2]."\t".$line[3]."\n";
			}
			$links{$rscaff}{$qscaff}++;
			$rlinks{$qscaff}{$rscaff}++;
			$links_weight{$rscaff}{$qscaff} += abs($line[0] - $line[1]);
			$rlinks_weight{$qscaff}{$rscaff} += abs($line[2] - $line[3]);
			
			my $this_link_swap_weight = ( ($line[0] - $line[1]) * ($line[2] - $line[3]) );
			if ($flipquery) { $links_swap_weight{$qscaff} += $this_link_swap_weight; }
			if ($flipreference)	{ $links_swap_weight{$rscaff} += $this_link_swap_weight; }
			$link_count++;
		}
	}

}
close CO;
close LI;
if ($debug) { 
	print STDERR "written a total of ".$link_count." links\n";
	print STDERR "stored info on ".(scalar keys %contigs2len_g1)." reference contigs\n";
	print STDERR "stored info on ".(scalar keys %contigs2len_g2)." query contigs\n";
	print STDERR "links counts are:\n";
	foreach my $k (keys %links) { 
		foreach my $l (keys %{$links{$k}}) { 
			print STDERR "\t".$k."->".$l." ".$links{$k}{$l}."\n";
		}	
	}
}
# foreach my $c (keys %contigs2len_g1) { 
	# print STDERR $contigs2len_g1{$c}."\t".$c."\n";
# }

#foreach my $d (keys %contigs2len_g2) { 
#	print STDERR $contigs2len_g2{$d}."\t".$d."\n";
#}

if ($debug && ( $flipreference | $flipquery))  {
	print STDERR "link 'flip weights' are:\n";
	foreach my $i (keys %links_swap_weight) { 
		print STDERR $i."\t".$links_swap_weight{$i}."\n";
	}
}
###############
#
# CHOOSE SOME COLORS EQUALLY SPACED AROUND HSV COLOUR SPACE. ALLOWS
# LINKS TO BE COLOURED WITH A SPECTRUM OF S,V FOR A PARTICULAR H, FOR SOME NICE EFFECTS..
#
###########
#by default this is an array of arrays, with # of colours as first index and then r,g,b values
#if gradual_shading > 1 then this is a an array of array of arrays, with # of intensities (the value of 'gradual shading') as first index,
#then number of colours,  then r,g,b
my @colors;
if ($gradual_shading > 1) { 
	
	my $gap_s = ($col_end_s - $col_start_s) / ($gradual_shading - 1);
	my $s_here = $col_start_s;
	my $gap_v = ($col_end_v - $col_start_v) / ($gradual_shading - 1);
	my $v_here = $col_start_v;
	foreach (my $i = 0; $i != $gradual_shading; $i++ ) {
		my @colorslice = ();
   	    	if ($colour_links_by_query) {
                	 @colorslice = colourset( ( scalar keys %contigs2len_g2 ) + 1 ,'equal_spacing',$col_start_h,$s_here,$v_here);
        	} else {
        	         @colorslice = colourset( ( scalar keys %contigs2len_g1 ) + 1,'equal_spacing',$col_start_h,$s_here,$v_here);
        	}
		push(@colors,\@colorslice);
		$s_here += $gap_s;
		$v_here += $gap_v;
	}
	if ($debug) { print STDERR "gradual_shading: generated a color array of ".(scalar @colors)." intensities and ".(scalar @{$colors[0]})." shades\n"; }
} else { 
	if ($colour_links_by_query) { 			
		@colors = colourset( ( scalar keys %contigs2len_g2 ) + 1 ,'equal_spacing',$col_start_h,$col_start_s,$col_start_v);
	} elsif ($colour_rainbow_links > 0 ) { 
		@colors = colourset( ( $colour_rainbow_links) ,'equal_spacing',$col_start_h,$col_start_s,$col_start_v);
	} else {
		@colors = colourset( ( scalar keys %contigs2len_g1 ) + 1,'equal_spacing',$col_start_h,$col_start_s,$col_start_v);
	}
	if ($debug) { print STDERR "no shading: generated a color array of ".(scalar @colors)." shades\n"; }
}

############
#ALL OF THE JUNK BELOW ALLOWS ORDERNG OF QUERY AND REFERENCE CONTIGS BY THE ORDER IN A TEXT OR FASTA FILE
############
#foreach my $c (@colors) { 
#	print join(" ",@$c)."\n";
#}

sub is_fasta_file { 
	my $fn = lc(shift);
	if ($fn =~ m/\.fa$/ || $fn =~ m/\.fasta$/ || $fn =~ m/\.fas$/ ) { return 1; } else { return 0; } 
}

#specially for nucmer, lets chop the rest of the names off
sub read_fasta_file { 
	my $fn = shift;
	my $title = "";
	my $seq = "";
	my %seqs;
	tie(%seqs, 'Tie::IxHash');
	open(F,"<$fn") or die("cannot read from fasta file ".$fn);
	while (<F>) { 
		chomp;	
		if (m/^>/) { 
			if ($seq ne "") { 
				$seqs{$title} = $seq;
			} 
			$seq = "";
			($title = $_) =~ s/^>//;
			$title =~ s/\s+.*$//g;
		} else { $seq .= $_; }
	}
	close F;
	$seqs{$title} = $seq;
	return \%seqs;
}

sub read_text_file { 
	my $fn = shift;
        my %lines;
	tie (%lines, 'Tie::IxHash');
	open(F,"<$fn") or die("cannot read from text file ".$fn);
	while (<F>) { 
		chomp;
		$_ =~ s/\s+.*$//g;
		$lines{$_}++;
	}
	close F;
	return \%lines;
}

sub order_by_hash_titles {
	my $local_debug = 0;
	my $list_to_order = shift;
	
	my @list_to_order = @{$list_to_order};
	my $ordering_hash = shift;
	my %ordering_hash;
	tie(%ordering_hash,'Tie::IxHash');
	%ordering_hash = %{$ordering_hash};	
	if ($local_debug) { print STDERR "sorting_hash keys are ".join(",",keys %ordering_hash)."\n"; }
	foreach my $l (@list_to_order) { if (! exists $ordering_hash{$l} ) { die("cannot find entry ".$l." in sort-against-hash\n"); } }
	my %ordering_key;
	my $i = 0;
	foreach my $k (keys %ordering_hash) { 
		$ordering_key{$k} = $i++;
	}
	if ($local_debug) { 
		print STDERR "sorting ".join(",",@list_to_order)." along:\n";
		foreach my $k (keys %ordering_key) { 
			print STDERR "\t".$k.":".$ordering_key{$k}."\n";
		}
	}
	my @out_list = sort { $ordering_key{$a} <=> $ordering_key{$b} } @list_to_order;
	if ($local_debug) { print STDERR "sorted ".join(",",@out_list)."\n";}
	return @out_list;	
}

my @contigs_g1_ordered = keys %contigs2len_g1;
if (! $ref_order ) { 
	@contigs_g1_ordered = sort { lc($a) cmp lc($b) } @contigs_g1_ordered;
} elsif (lc($ref_order) ne "rel" && lc($ref_order) ne "relw") { 
	#ref_order is presumably a filename
     	my $seqs;
	if (is_fasta_file($ref_order) ) {
                if ($debug) { print "reading order for ref contigs from fasta file ".$ref_order."\n"; }
       		$seqs = read_fasta_file($ref_order);
		if ($debug) { print " read ".(scalar keys %$seqs)." unique seq titles from fasta file\n"; } 
	} else { 
		if ($debug) { print "reading order for ref contigs from text file ".$ref_order."\n"; }
                $seqs = read_text_file($ref_order);
                if ($debug) { print " read ".(scalar keys %$seqs)." unique seq titles from text file\n"; }
	}
	if ($parse_contig_names) { 
		my %new_seqs;
		tie(%new_seqs,'Tie::IxHash');
		foreach my $k (keys %$seqs) {
#			if ($debug) { print STDERR "looking to insert ".$k." in seqs\n"; }
			if (exists  $r_scaff_translation{$k} ) { 
				$new_seqs{$r_scaff_translation{$k}}++;
#				if ($debug) { print STDERR " inserting ".$r_scaff_translation{$k}." instead\n"; }
			}
#			} else { if ($debug) { print STDERR "not found \n"; } }
		}
		$seqs = \%new_seqs;
	}
	if ($debug) { print STDERR "sort_against_hash is ".join(",",keys %$seqs)."\n"; } 
	@contigs_g1_ordered = order_by_hash_titles(\@contigs_g1_ordered,$seqs);
	if ($debug) { print STDERR "final sorted order is ".join(",",@contigs_g1_ordered)."\n"; } 
}
my @contigs_g2_ordered = keys %contigs2len_g2;
if (! $query_order ) { 
	@contigs_g2_ordered = sort { lc($a) cmp lc($b) } @contigs_g2_ordered;
} elsif (lc($query_order) ne "rel" && lc($query_order) ne "relw") { 
	#query_order is presumably a filename
	my $seqs;
	if (is_fasta_file($query_order) ) { 
		if ($debug) { print "reading order for query contigs from fasta file ".$query_order."\n"; } 
		$seqs = read_fasta_file($query_order);
                if ($debug) { print " read ".(scalar keys %$seqs)." unique seq titles from fasta file\n"; }
	} else { 
		if ($debug) { print "reading order for query contigs from text file ".$query_order."\n"; }
                $seqs = read_text_file($query_order);
                if ($debug) { print " read ".(scalar keys %$seqs)." unique seq titles from text file\n"; }
	}
 	if ($parse_contig_names) {
                my %new_seqs;
                tie(%new_seqs,'Tie::IxHash');
                foreach my $k (keys %$seqs) {
                        if (exists  $q_scaff_translation{$k} ) {
                                $new_seqs{$q_scaff_translation{$k}}++;
                        }
                }
                $seqs = \%new_seqs;
        }
 	if ($debug) { print STDERR "sort_against_hash is ".join(",",keys %$seqs)."\n"; }
	@contigs_g2_ordered = order_by_hash_titles(\@contigs_g2_ordered,$seqs);
	@contigs_g2_ordered = reverse @contigs_g2_ordered;
	if ($debug) { print STDERR "final sorted order is ".join(",",@contigs_g2_ordered)."\n"; }
	
}

if (lc($ref_order) eq "rel" || lc($ref_order) eq "relw") { 
	my %q_order_now;
	my $i = 0;
	foreach my $k (@contigs_g2_ordered) { 
		$q_order_now{$k} = $i++;
	} 
	my %total_links;
	my %total_weight;
	foreach my $l (@contigs_g1_ordered) { 
		
		foreach my $k (keys %{$links{$l}} ) { 
			if ( lc($ref_order) eq "rel" ) { 
				$total_links{$l} += $links{$l}{$k};
				$total_weight{$l} += $links{$l}{$k} * $q_order_now{$k};
			} 
			if ( lc($ref_order) eq "relw") { 
				$total_links{$l} += $links_weight{$l}{$k};
				$total_weight{$l} += $links_weight{$l}{$k} * $q_order_now{$k};
			}
		}
		if ($debug) { print STDERR "\t".$l." totalw=".$total_weight{$l}." totall=".$total_links{$l}."\n"; }
	}
	foreach my $k (keys %total_weight) { $total_weight{$k} = $total_weight{$k} / $total_links{$k}; }
	if ($debug) { 
		print "order here is ".join(",",@contigs_g1_ordered)."\n";
	}
	@contigs_g1_ordered = sort { $total_weight{$a} <=> $total_weight{$b} } @contigs_g1_ordered;
	@contigs_g1_ordered = reverse @contigs_g1_ordered;
	if ($debug) { print "order now is ".join(",",@contigs_g1_ordered)."\n"; } 
}
if (lc($query_order) eq "rel" || lc($query_order) eq "relw") { 
        my %q_order_now;
        my $i = 0;
        foreach my $k (@contigs_g1_ordered) {
                $q_order_now{$k} = $i++;
        }
        my %total_links;
        my %total_weight;
        foreach my $l (@contigs_g2_ordered) {
                
                foreach my $k (keys %{$rlinks{$l}} ) {  
					if ( lc($query_order) eq "rel" ) {
							$total_links{$l} += $rlinks{$l}{$k};
							$total_weight{$l} += $rlinks{$l}{$k} * $q_order_now{$k};
					}
					if ( lc($query_order) eq "relw" ) {
							$total_links{$l} += $rlinks_weight{$l}{$k};
							$total_weight{$l} += $rlinks_weight{$l}{$k} * $q_order_now{$k};
					}
                }
                if ($debug) { print STDERR "\t".$l." totalw=".$total_weight{$l}." totall=".$total_links{$l}."\n"; }
        }
        foreach my $k (keys %total_weight) { $total_weight{$k} = $total_weight{$k} / $total_links{$k}; }
        if ($debug) { 
                print "order here is ".join(",",@contigs_g2_ordered)."\n";
        }
        @contigs_g2_ordered = sort { $total_weight{$a} <=> $total_weight{$b} } @contigs_g2_ordered;
	@contigs_g2_ordered = reverse @contigs_g2_ordered;
        if ($debug) { print "order now is ".join(",",@contigs_g1_ordered)."\n"; } 

} 


###########
#
# DONE SORTING - NOW PRINT KARYOTYPE FILE IN APPROPRIATE ORDER, AND WITH APPROPRIATE COLORS
#
###########

open(DATA,">",$prefix."karyotype.txt");
my $i = 1;

if ($debug) { 
	print STDERR "now have ".(scalar @contigs_g1_ordered)." ordered ref contigs and ".(scalar @contigs_g2_ordered)." ordered query contigs\n";
	print STDERR "ref order is ".join(",",@contigs_g1_ordered).", query order is ".join(",",@contigs_g2_ordered)."\n"; 
}
foreach my $k (@contigs_g1_ordered) { 
#foreach my $k (keys %contigs2len_g1) { 
	if ($colour_links_by_query) { 
		if ($no_ref_labels) { 
			print DATA "chr - ".$k." £no-print£ 0 ".$contigs2len_g1{$k}." c0";
		} 
		elsif ($parse_contig_names) {
			print DATA "chr - ".$k." ".$r_scaff_translation_back{$k}." 0 ".$contigs2len_g1{$k}." c0";	
		} else {
			print DATA "chr - ".$k." ".$k." 0 ".$contigs2len_g1{$k}." c0";
		}
		print DATA "\n";
	#	print DATA "chr - ".$k." ".$i++." 0 ".$contigs2len_g1{$k}." c0\n";
		#$i++;
	} else { 
		if ($no_ref_labels) { 
			print DATA "chr - ".$k." £no-print£ 0 ".$contigs2len_g1{$k}." c".$i;
		}
		elsif ($parse_contig_names) {
			print DATA "chr - ".$k." ".$r_scaff_translation_back{$k}." 0 ".$contigs2len_g1{$k}." c".$i;
		} else { 
			print DATA "chr - ".$k." ".$k." 0 ".$contigs2len_g1{$k}." c".$i;
			
		}
	#	print DATA "chr - ".$k." ".$i++." 0 ".$contigs2len_g1{$k}." c".$i."\n";
		$i++;
        	if ($gradual_shading > 1) {
        	        print DATA ".".(ceil($gradual_shading / 2))."\n";
        	} else {
        	        print DATA "\n";
	        }
	}
}
my $j = 1;
foreach my $k1 (@contigs_g2_ordered) { 
#foreach my $k1 (keys %contigs2len_g2) { 
	if ($colour_links_by_query) { 
		if ($no_query_labels) { 
			print DATA "chr - ".$k1." £no-print£ 0 ".$contigs2len_g2{$k1}." c".$j;
		}
		elsif ($parse_contig_names) {
                	print DATA "chr - ".$k1." ".$q_scaff_translation_back{$k1}." 0 ".$contigs2len_g2{$k1}." c".$j;
                } else {
			print DATA "chr - ".$k1." ".$k1." 0 ".$contigs2len_g2{$k1}." c".$j;
		}
#		print DATA "chr - ".$k1." ".$j++." 0 ".$contigs2len_g2{$k1}." c".$j."\n";
		$j++;
        	if ($gradual_shading > 1) {
        	        print DATA ".".(ceil($gradual_shading / 2))."\n";
        	} else {
        	        print DATA "\n";
        	}
	} else {
		if ($no_query_labels) { 
			print DATA "chr - ".$k1." £no-print£ 0 ".$contigs2len_g2{$k1}." c0";
		}
		elsif ($parse_contig_names) {
                
			print DATA "chr - ".$k1." ".$q_scaff_translation_back{$k1}." 0 ".$contigs2len_g2{$k1}." c0";
		} else { 
			print DATA "chr - ".$k1." ".$k1." 0 ".$contigs2len_g2{$k1}." c0";
		} 
#		print DATA "chr - ".$k1." ".$j++." 0 ".$contigs2len_g2{$k1}." c0\n";
		#$j++;
		print DATA "\n";
	}
}
close DATA;

####
#
# NOW PRINT MAIN circos.conf file.
# this includes colour definitions and 'rules' to define colours of links, 
# so is quite long
#
###

open(CONF,">circos.conf");
print CONF "karyotype = ".$prefix."karyotype.txt\n";
print CONF "\n";
if ($flipquery || $flipreference) { 
	my @flippers;
	my @all;
	#if ($flipreference) { 
		foreach my $i ( @contigs_g1_ordered) { 
			push(@all,$i);
		}
#	}
	#if ($flipquery) { 
		foreach my $i ( @contigs_g2_ordered) { 
			push(@all,$i);
		}
#	}
	foreach my $i (keys %links_swap_weight) { 
		if ($links_swap_weight{$i} > 0) { push(@flippers,$i); }
	}
	print CONF "chromosomes = ".join(";",@all)."\n";
	print CONF "chromosomes_reverse = ".join(";",@flippers)."\n";
}
print CONF "\n";
print CONF "<ideogram>\n";
print CONF "\n";
print CONF "<spacing>\n";
print CONF "default = 0.005r\n";
print CONF "</spacing>\n";
print CONF "\n";
print CONF "radius    = 0.9r\n";
print CONF "thickness = 20p\n";
print CONF "fill      = yes\n";
print CONF "show_label = yes\n";
print CONF "label_with_tag = no\n";
print CONF "label_font = light\n";
print CONF "label_radius = 1r + 4p\n";
print CONF "label_center = yes\n";
print CONF "label_size     = ".$label_size_param."\n";
#print CONF "label_color    = grey\n";
print CONF "label_parallel = yes\n";
print CONF "label_case     = upper \n";
#print CONF "label_format   = eval(sprintf(\"chr%s\",var(label)))\n";
#change '£no-print£' into blank
print CONF "label_format = eval(my \$x = var(label); \$x eq \"£no-print£\" ? \"\" : \$x)\n";
print CONF "\n";


print CONF "</ideogram>\n";
print CONF "\n";
if ($zoom_query) { 
	my @magnificat;
	foreach my $k1 (@contigs_g2_ordered) { 
		push(@magnificat,$k1.":".$zoom_query);
	}
	print CONF "chromosomes_scale = ".join(";",@magnificat)."\n";
#	print CONF "<zooms>\n";
#	foreach my $k1 (@contigs_g2_ordered) {
#		print CONF "<zoom>\n";
#		print CONF "chr = ".$k1."\n";
#		print CONF "start = 0\n";
#		print CONF "end = ".$contigs2len_g2{$k1}."\n";
#		print CONF "scale = ".$zoom_query."\n";
#		print CONF "</zoom>\n";
#	} 
#	print CONF "</zooms>\n";

}
print CONF "<colors>\n";
if ($gradual_shading > 1) { 
	for (my $i = 0; $i != scalar @colors; $i++) { 
		for (my $j = 0; $j != scalar @{$colors[$i]}; $j++) { 
			if ($j == 0) { 
				#only need one 'intensity' for c0
				if ($i == ceil($gradual_shading / 2) ) { 
					print CONF "c".$j." = ".join(",",@{$colors[$i][$j]})."\n";
				} 
			} else { 
				print CONF "c".$j.".".$i." = ".join(",",@{$colors[$i][$j]})."\n";
			}
		}
	}
} else  {
	for (my $i = 0; $i != scalar @colors; $i++) { 
		print CONF "c".$i." = ".join(",",@{$colors[$i]})."\n";
	}
}
print CONF "</colors>\n";
print CONF "<links>\n";
print CONF "<link>\n";
if ($ribbons) { 
	print CONF "ribbon       = yes\n";
}
print CONF "file          = ".$prefix."links.txt\n";
#print CONF "color         = black_a5\n";
print CONF "radius        = 0.95r\n";
print CONF "bezier_radius = 0.1r\n";
print CONF "thickness     = 1\n";

#generate 'cut points' for intensity shading on each of the 
#query / ref contigs
#


print CONF "<rules>\n";
if ($colour_links_by_query) { 
	my $i = 1;
	foreach my $k (@contigs_g2_ordered) { 
	#foreach my $k (keys %contigs2len_g2) { 
		if ($gradual_shading > 1) {
			my $len = $contigs2len_g2{$k};
			my $cut_start = 0;
			my $cut_sizes = $len / $gradual_shading;
			if($debug) { 
				my @cut_sites = ($cut_start);
				foreach (my $i =0; $i != $gradual_shading; $i++) { 
					push(@cut_sites,floor($cut_start += $cut_sizes));
				}
		  		 print STDERR "for contig ".$k." (length ".$len.") shading cuts are at ".join(",",@cut_sites)."\n"; 
				$cut_start = 0;
			}
		       
			foreach (my $j = 0; $j != $gradual_shading; $j++) { 
  				print CONF "<rule>\n";
                        	print CONF "condition = var(chr2) eq \"".$k."\"\n";
				print CONF "condition = on(".$k.",".floor($cut_start).",".floor($cut_start+$cut_sizes).")\n";
				$cut_start += $cut_sizes;
				print CONF "color = c".$i.".".$j."\n";
        	        	print CONF "z = ".$i."\n";
                		print CONF "</rule>\n";			
			}
			$i++;
			#print CONF "color = c".$i++.".".(ceil($gradual_shading / 2))."\n";
        	} else {
			print CONF "<rule>\n";
                	print CONF "condition = var(chr2) eq \"".$k."\"\n";
        	        print CONF "color = c".$i++."\n";
			print CONF "z = ".$i."\n";
                	print CONF "</rule>\n";
		}
		##print CONF "color = c".$i++."\n";
	}
} elsif ($colour_rainbow_links > 0) { 
	print CONF "<rule>\n";
	# always trigger this rule
	print CONF "condition  = 1\n";
	# use the link's value to sample from a list of colors
	print CONF "color      = eval((qw(";
	my @col_list;
	foreach (my $i =0; $i != $colour_rainbow_links; $i++) { 
		push(@col_list,"c".$i);
	}
	print CONF join(" ",@col_list)."))[ var(value) ])\n";
	# continue parsing other rules
	print CONF "flow       = continue\n";
	print CONF "</rule>\n";

} else { 
        my $i = 1;
       foreach my $k (@contigs_g1_ordered) { 
	# foreach my $k (keys %contigs2len_g1) {
		if ($gradual_shading > 1) { 
     			my $len = $contigs2len_g1{$k};
                        my $cut_start = 0;
                        my $cut_sizes = $len / $gradual_shading;
                        if($debug) {
                                my @cut_sites = ($cut_start);
                                foreach (my $i =0; $i != $gradual_shading; $i++) {
                                        push(@cut_sites,floor($cut_start += $cut_sizes));
                                }
                                 print STDERR "for contig ".$k." (length ".$len.") shading cuts are at ".join(",",@cut_sites)."\n";
                                $cut_start = 0;
                        }

                        foreach (my $j = 0; $j != $gradual_shading; $j++) {
                                print CONF "<rule>\n";
                               	print CONF "condition = var(chr1) eq \"".$k."\"\n";
                                print CONF "condition = on(".$k.",".floor($cut_start).",".floor($cut_start+$cut_sizes).")\n";
                                $cut_start += $cut_sizes;
                                print CONF "color = c".$i.".".$j."\n";
                                print CONF "z = ".$i."\n";
                                print CONF "</rule>\n";
                        }
                        $i++;
		} else { 

			print CONF "<rule>\n";
	                print CONF "condition = var(chr1) eq \"".$k."\"\n";
	                print CONF "color = c".$i++."\n";
			print CONF "z = ".$i."\n";
			print CONF "</rule>\n";
        	}
	}
}
#print CONF "<rule>\n";
#print CONF "condition = var(interchr)\n";
#print CONF "color = red\n";
#print CONF "</rule>\n";
print CONF "</rules>\n";
print CONF "</link>\n";
print CONF "</links>\n";


print CONF "<image>\n";
print CONF "file*  = ".$prefix."circos.png\n";
print CONF "\n";
print CONF "# Included from Circos distribution.\n";
print CONF "<<include etc/image.conf>>\n";
print CONF "</image>\n";
print CONF "\n";
print CONF "# RGB/HSV color definitions, color lists, location of fonts, fill patterns.\n";
print CONF "# Included from Circos distribution.\n";
print CONF "<<include etc/colors_fonts_patterns.conf>>\n";
print CONF "\n";
print CONF "# Debugging, I/O an dother system parameters\n";
print CONF "# Included from Circos distribution.\n";
print CONF "<<include etc/housekeeping.conf>>\n";
print CONF "\n";



######
#
# STOLEN CODE TO BISECT HSV COLOURSPACE, AND CONVERT TO RGB.
#
#####

=item colourset($num_colours, $method)
Function to grab a set of well spaced colours from an HSV colour wheel.
    $num_colours - The number of colour values to produce, must be greater than 0 but no bigger than 360
    $method - The method for selecting colours over HSV colour space, either 'equal_spacing' or for around 10 colours 'chroma_bisection' is better.
either pass 2 parameters, or specify H,S,V starting point.
Returns an array of RGB values of the form ([R,G,B]) and undef on $num_colours out of bounds
=cut


sub colourset {
        if (scalar @_ != 5) { 
		push(@_,0.0,0.65,1.0);
	}  
	my ($num_colours,$method,$start_H,$start_S,$start_V) = @_;
#should be 360
        if ($num_colours <= 0 or $num_colours > 720) {
                warn "Number of colours requested out of bounds.";
                return undef;
        }
        $method = 'chroma_bisection' unless $method;
        #Internal sub to randomly shuffle an array
        sub fisher_yates_shuffle {
                my ($array) = @_;
                my $current;
                for ($current = @$array; --$current; ) {
                        my $selected = int rand ($current+1);
                        next if $current == $selected;
                        #Reverse the slice between current position and the randomly selected
                        @$array[$current,$selected] = @$array[$selected,$current];
                }
                return $array;
        }
        
        #Colours to return
        my %colours;
        
        #Default Hue Saturation and Value, saturation of 0.65 gives a more pastel feel!
        my ($Hue, $Saturation, $Value) = ($start_H,$start_S,$start_V);
        
        #The interval to space colours around the wheel if equal
        my $hsv_interval = 360 / $num_colours;
        
        #Array of degrees for reuse to create ranged arrays with a given interval
        my @degrees = 1..360;
        
        #Iteratively bisect each chroma segment so that the first 6 colours are well spaced perceptually.
        #However after 12 colours we will have increasing pairs that are more confused as 
        #they are increasingly close to each other compared to the rest of the colours!
        #To get around this problem of selecting closely around a single bisection, we jump around the 
        #chroma randomly sampling.
        if ($method eq 'chroma_bisection') {
                #The current cycle of chroma bisection
                my $hsv_cycle = 1;
                #Number of colours selected by bisecting chroma so far
                my $colours_selected = 0;
                
                #While we still have colours to select
                while ($colours_selected != $num_colours) {
                        #Work out the size of interval to use this cycle around the wheel
                        $hsv_interval = 60 / $hsv_cycle;
                        #Get all the Hues for this cycle that haven't already been examined and are on the line of bisection
                        my @Hues = grep { (not $_ % $hsv_interval) && (not exists $colours{$_%360}) } @degrees;
                        #Shuffle so that we don't take from around the same chroma all the time, only perceptually worthwhile after 12th colour
                        fisher_yates_shuffle(\@Hues) if $hsv_cycle > 2;
                        
                        #While we still have hues to select from in this cycle
                        while (@Hues) {
                                #Finish if we have enough colours
                                last if $colours_selected == $num_colours;
                                #Consume a Hue from this cycle
                                $Hue = shift @Hues;
                                #360 should be 0 for red
                                $Hue %= 360;
                                #Store this Hue and mark selection
                                $colours{$Hue} = hsv2rgb($Hue,$Saturation,$Value) ;
                                $colours_selected++;
                        }
                        $hsv_cycle++;
                }
        }
        
        #Just space colours even distances apart over the HSV colour wheel.
        #You have slightly odd/garish colours coming out, but you dont get uneven perceptual distance
        #between pairs of colours. This scales far better despite the horrible colours.
        elsif ($method eq 'equal_spacing') {    
                foreach $Hue (1..$num_colours) {
                        $Hue = ($Hue * $hsv_interval) % 360;
                        $colours{$Hue} = hsv2rgb($Hue,$Saturation,$Value) ;
                }
        }
        
        #Otherwise return nothing and warn the programmer
        else {
                warn "Colourset method not known, use either 'equal_spacing' or for fewer colours 'chroma_bisection'";
                return undef;
        }
        
        #Shuffle final colours so that even if we do use chroma_bisection closer colours will hopefully not be sequential
        @_ = values %colours;
#        fisher_yates_shuffle(\@_);
        return @_;
}


=item hsv2rgb($Hue, $Saturation, $Value)
Function to convert HSV colour space values to RGB colour space.
Returns RGB value as [R,G,B]
=cut
sub hsv2rgb {
        my ($Hue,$Saturation,$Value) = @_;
        my ($Red,$Green,$Blue) = (0,0,0);
        
        #Check the input and warn if it's a bit wrong
        warn "Invalid Hue component of HSV colour passed, with value: $Hue." unless ($Hue >= 0.0 and $Hue <= 360.0);
        warn "Invalid Saturation component of HSV colour passed, with value: $Saturation." unless($Saturation >= 0.0 and $Saturation <= 1.0);
        warn "Invalid Value component of HSV colour passed, with value: $Value." unless ($Value >= 0.0 and $Value <= 1.0);
        
        #If colour has no saturation return greyscale RGB
        if ($Saturation == 0) {
                $Red = $Green = $Blue = $Value;
                return [$Red, $Green, $Blue];
        }
        
        #Partition the Hue into the 5 different colour chroma and then map each of these to RGB based on the colour theory
        $Hue /= 60.0;
        my $Chroma = floor($Hue) % 6; 
        my $H_d = $Hue - $Chroma; 
        
        #RGB cube components
        my ($I,$J,$K) = ( $Value * ( 1 - $Saturation ),
                                   $Value * ( 1 - $Saturation * $H_d ),
                                   $Value * ( 1 - $Saturation * ( 1 - $H_d ) )
                                    );
        
        #Map components to RGB values per chroma
        if ($Chroma == 0) { ($Red,$Green,$Blue) = ($Value,$K,$I); }
        elsif ($Chroma == 1) { ($Red,$Green,$Blue) = ($J,$Value,$I); }
        elsif ($Chroma == 2) { ($Red,$Green,$Blue) = ($I,$Value,$K); }
        elsif ($Chroma == 3) { ($Red,$Green,$Blue) = ($I,$J,$Value); }
        elsif ($Chroma == 4) { ($Red,$Green,$Blue) = ($K,$I,$Value); }
        else{ ($Red,$Green,$Blue) = ($Value,$I,$J); }
        
        #Return the RGB value in the integer range [0,255] rather than real [0,1]
        return [floor($Red * 255),floor($Green * 255),floor($Blue * 255)];
}
