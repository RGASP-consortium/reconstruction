2013, July 4th
RGASP_analyzer
(developed for RGASP, for details see: http://www.gencodegenes.org/rgasp/)

This small program can be used to compare an automatically generated gtf file with a
gtf file annotation. It will output how many features (exons, transcript and genes) are
within the annotation and the prediction and how many of the predicted features are correct. 

Installation:
	- make sure Java 1.5 or newer is installed
	- jargs.jar must be on the Java classpath
	- bash must be available
	- awk must be installed and in the path

Usage: java -Xms1G -Xmx2G org.rgasp.PredictionAnalyser
	-l, --level			... one of "cds", "exon", "both"
	
							if "cds" is specified, only cds records within the annotation and
							the prediction are considered
							
							if "exon" is specified only exon records within the annotation and
							the prediction are considered
							
							if "both" is specified a "cds" and an "exon" analysis
							will be performed separately
							
	-m, --mode			... one of "fixed", "flexible", "both"
	
							if "fixed" mode is specified, predicted records are verified
							only if there is an annotation record with the same borders
							
							if "flexible" mode is specified, predicted records, that correspond to
							a first or last exon of a transcript in the annotation, may have flexible
							outer borders
							
							if "both" is specified, fixed and flexible mode analysis will be performed
							separately
							
	-a, --annotation	... path to annotation file
	-p, --prediction	... path to prediction file / folder
	-o, -- output		... output folder
	
Example:
	java -Xms1G -Xmx2G -jar RGASP_analyzer.jar -l exon -m both -a path2annotation -p path2prediction -o outputfolder

File format:
	Both the annotation and the prediction files are expected to be in gtf file format.
	Please note that the key-value pairs gene_id and transcript_id are required!
	(for details see: http://www.sanger.ac.uk/resources/databases/encode/gencodeformat.html)
	
	Please note, that this analysis only considers gtf records with the feature-type exon or CDS.
	Make sure, that you adapt your gtf file accordingly.
	Example: If your gtf file contains information about coding exons and UTRs, use this information
	to create exon records.
	
	Please make sure, that the file format for the annotation and the prediction files are consistent,
	especially make sure that:
		- chromosome names are called consistent
		- stop codons are consistently within the coding exons or not

Output:
	In the folder specified with the -o, --output parameter the following tables will be created
	
	- LEVEL_MODE_statistics.txt
	  lists for a single prediction or each prediction file in the prediction folder how many exons,
	  transcripts, or genes have been predicted, and how many of them are true positives.
	  Example:
	  	sensitivity for (coding) exons can be calculated by: TP_Exons / Num_Exons_Ref
	  	specificity for (coding) exons can be calculated by: TP_Exons_Pred / Num_Exons_Pred
	
	- LEVEL_MODE_exon.tbl
	  binary table that lists for each exon of an annotated transcript whether it has been
	  predicted (1) or not (0) by a prediction file.
	  Exons are identified by their gene and transcript identifier and their position within the transcript
	
	- LEVEL_MODE_exon_unique.tbl
	  binary table as LEVEL_MODE_exon.tbl except that exons that are shared between different transcript
	  isoforms will only appear once
	
	- LEVEL_MODE_transcript.tbl
	  binary table that lists for each transcript of an annotation file whether it has been predicted (1) or not (0)
	  by a prediction file.
	
	- LEVEL_MODE_gene.tbl
	  binary table that lists for each gene of an annotation file whether it has been predicted (1) or not (0)
	  by a prediction file.
	  
	While the statistic table gives you fast access to sensitivity and specificity values, the other tables can
	be used to make more complex follow-up analyses.
	

Memory: 
	The memory Java allocated by default will most probably not be sufficient to run
	your analysis. Therefore, you would usually run the analysis with the -Xms1G -Xmx2G arguments.
	Depending on your input file sizes, you might want to increase these argument.
	
Contact: When facing any problems please contact:
	Tamara Steijger
	European Bioinformatics Institute
	tamara.steijger@ebi.ac.uk
