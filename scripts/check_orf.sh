#!/bin/bash

# This script splits a CDS sequence in FASTA format into codons, highlights in-frame stop codons,
# and indicates "ORF COMPLETE" if the sequence has no in-frame stop codons.
# $1 should be provided as the FASTA file, for example: 'fasta.fa'.

fa=$1

# Process the multi-FASTA file to handle headers and concatenate sequence lines.
awk '/^>/{if(NR>1)print ""; print; next} {printf "%s", $0} END{print ""}' $fa |

# Remove hyphens from sequences.
sed 's/-//g' |

# Remove the last three characters from sequence lines.
sed '/>/! s/.\{3\}$//' |

# Insert tabs every three characters in sequence lines.
sed '/>/! s/.\{3\}/&\t/g' |

# Search for stop codons "taa", "tag", or "tga" (case-insensitive).
# If no stop codons are found, output "ORF COMPLETE".
egrep -i --color=auto "taa|tag|tga" || echo "ORF COMPLETE"
