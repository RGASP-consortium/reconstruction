package org.rgasp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.rgasp.GTFRecord.Strand;
import org.rgasp.PredictionAnalyser.Feature;

public class GFFParser {
	
	protected BufferedReader file_reader;
	protected long m_size;
	protected long m_currentOffset;
	protected Feature feature;
	
	public GFFParser (String filename, Feature feature) throws IOException{
		System.out.println("Parsing " + filename);
		File temp = new File(filename);
		m_size = temp.length();
		
		FileReader file = new FileReader(filename);
		file_reader = new BufferedReader(file);
		this.feature = feature;
	}
	
	private String[] split(String line){
		String[] temptokens = line.split("[;\t\\s\"]");
		ArrayList<String> tokens = new ArrayList<String>();
		//remove empty entries
		for(int i = 0; i < temptokens.length; ++i) {
			  if(temptokens[i] != null && !temptokens[i].equals("")) {
				  tokens.add(temptokens[i]);
			  }
		}
		String[] copytokens = new String[tokens.size()];
		return (String[]) (tokens.toArray(copytokens));
	}

	public void getTranscriptList(Map<String, List<GTFRecord>> exon_list, Map<String, Map<String, GTFTranscript>> transcript_list,
			Map<String, Map<String, GTFGene>> gene_list) throws IOException {
		Map<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>> unique_exons = new HashMap<String, Map<Integer, Map<Integer, Map<Strand, GTFRecord>>>>();
		
		String line = "";
		while (file_reader.ready()) {
			line = file_reader.readLine();
			if (line.equals("") || line.startsWith("#")) continue;
			String[] tokens = split(line);
			if (tokens == null || tokens.length < 7) System.out.println("Error! Missing columns in line:\n" + line);
			String sequence = tokens[0];
			String type = tokens[2].toUpperCase();
			if (type.equals("CODING_EXON")) type = "CDS";
			if (type.startsWith("CDS")) type = "CDS";
			if (!(type.equals("EXON") || type.equals("CDS"))) continue;
			if (feature == Feature.CDS && !type.equals("CDS")) continue;
			if (feature == Feature.EXON && !type.equals("EXON")) continue;
			int start_pos = Integer.parseInt(tokens[3]);
			int end_pos = Integer.parseInt(tokens[4]);
			String strand;
			if (tokens[6].equals("+")) strand = "Strand_Positive";
			else if (tokens[6].equals("-")) strand = "Strand_Negative";
			else strand = "Strand_Unknown";
			String transcript_id = "";
			String gene_id = "";
			if (tokens[8].equals("transcript_id")) {
				transcript_id = tokens[9];
				gene_id = tokens[11];
			}
			else if (tokens[10].equals("transcript_id")) {
				transcript_id = tokens[11];
				gene_id = tokens[9];
			}
			
			transcript_id = gene_id + ":" + transcript_id;
			GTFRecord record = new GTFRecord(sequence, type, start_pos, end_pos, strand, transcript_id);
			GTFTranscript transcript = null;
			
			
			if (!unique_exons.containsKey(sequence)) unique_exons.put(sequence, new HashMap<Integer, Map<Integer, Map<Strand, GTFRecord>>>());
			if (!unique_exons.get(sequence).containsKey(start_pos)) unique_exons.get(sequence).put(start_pos, new HashMap<Integer, Map<Strand, GTFRecord>>());
			if (!unique_exons.get(sequence).get(start_pos).containsKey(end_pos)) unique_exons.get(sequence).get(start_pos).put(end_pos, new HashMap<Strand, GTFRecord>());
			if (!unique_exons.get(sequence).get(start_pos).get(end_pos).containsKey(Strand.valueOf(strand))) {
				unique_exons.get(sequence).get(start_pos).get(end_pos).put(Strand.valueOf(strand), record);
				if (!exon_list.containsKey(sequence)) exon_list.put(sequence, new ArrayList<GTFRecord>());
				exon_list.get(sequence).add(record);
			}
			else {
				record = unique_exons.get(sequence).get(start_pos).get(end_pos).get(Strand.valueOf(strand));
				record.addContainer(transcript_id);
			}
			if (!transcript_list.containsKey(sequence)) transcript_list.put(sequence, new HashMap<String, GTFTranscript>());
			if (!transcript_list.get(sequence).containsKey(transcript_id)) {
				transcript = new GTFTranscript(record, gene_id);
				GTFTranscript.existing_transcripts.put(transcript_id, transcript);
				transcript.transcript_id = transcript_id;
				transcript_list.get(sequence).put(transcript_id, transcript);
			}
			else {
				transcript_list.get(sequence).get(transcript_id).addExon(record);
			}
			if (gene_list != null) {
				if (!gene_list.containsKey(sequence)) gene_list.put(sequence, new HashMap<String, GTFGene>());
				if (!gene_list.get(sequence).containsKey(gene_id)) {
					GTFGene gene = new GTFGene(transcript);
					gene.gene_id = gene_id;
					gene_list.get(sequence).put(gene_id, gene);
				}
				else if (transcript != null) {
					gene_list.get(sequence).get(gene_id).addTranscript(transcript);
				}
			}
		}
		for (List<GTFRecord> list : exon_list.values()) {
			Collections.sort(list);
		}
		for (String sequence : transcript_list.keySet()) {
			for (GTFTranscript transcript : transcript_list.get(sequence).values()) {
				Collections.sort(transcript.getExons());
			}
		}
		file_reader.close();
		file_reader = null;
	}

}
