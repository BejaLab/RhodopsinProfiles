The rhodopsin sequences are from https://github.com/BejaLab/RhodopsinsReview, exluding 8TMs and HeRs:

	git clone https://github.com/BejaLab/RhodopsinsReview
	sed -f rename_header.sed < RhodopsinsReview/data/rhodopsins.fasta | seqkit grep -nrvf skip.txt > sequences.fasta
