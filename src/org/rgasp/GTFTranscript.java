package org.rgasp;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.rgasp.GTFRecord.Strand;
import org.rgasp.PredictionAnalyser.Level;

public class GTFTranscript extends GTFRecord {
	
	protected static Map<String, GTFTranscript> existing_transcripts = new HashMap<String, GTFTranscript>();
	protected String transcript_id;
	
	protected List<GTFRecord> exons;
	
	public List<String> isPredicted(Level level) {
		ArrayList<String> matching_pred_transcripts = new ArrayList<String>();
		//check single exon transcripts first
		if (level == Level.FLEXIBLE && exons.size() == 1) {
			if (exons.get(0).overlap < 0.9) return new ArrayList<String>();
			boolean single_exon_prediction = false;
			for (GTFRecord overlapping_exon : exons.get(0).overlapingExons) {
				for (String s : overlapping_exon.container) {
					if (existing_transcripts.get(s).getExons().size() == 1) {
						single_exon_prediction = true;
						matching_pred_transcripts.add(s);
					}
				}
			}
			if (single_exon_prediction) return(matching_pred_transcripts);
			else new ArrayList<String>();
		} else if (level == Level.FIXED && exons.size() == 1) {
			if (exons.get(0).matchingExons.size() < 1) return new ArrayList<String>();
			boolean single_exon_prediction = false;
			for (GTFRecord matching_exon : exons.get(0).matchingExons) {
				for (String s : matching_exon.container) {
					if (existing_transcripts.get(s).getExons().size() == 1) {
						single_exon_prediction = true;
						matching_pred_transcripts.add(s);
					}
				}
			}
			if (single_exon_prediction) return(matching_pred_transcripts);
			else new ArrayList<String>();
		}
		// for longer transcripts get candidate transcripts, e.g., check transcripts of the first exon (fixed or flexible)
		List<String> candidate_list;
		if (level == Level.FIXED && exons.get(0).matchingExons.size() > 0) candidate_list = exons.get(0).matchingExons.get(0).container;
		else {
			if (exons.get(0).matchingExons.size() > 0) candidate_list = exons.get(0).matchingExons.get(0).container;
			else candidate_list = new ArrayList<String>();
			if (PredictionAnalyser.end_positions.containsKey(exons.get(0).sequence) && PredictionAnalyser.end_positions.get(exons.get(0).sequence).containsKey(exons.get(0).strand)
					&& PredictionAnalyser.end_positions.get(exons.get(0).sequence).get(exons.get(0).strand).contains(exons.get(0).end_pos)) {
				for (GTFRecord rec : PredictionAnalyser.end_positions_map.get(exons.get(0).sequence).get(exons.get(0).strand).get(exons.get(0).end_pos)) {
					for (String s : rec.container) {
						if (!candidate_list.contains(s)) candidate_list.add(s);
					}
				}
			}
		}
		for (String candidate_transcript : candidate_list) {
			if (does_transcript_match(candidate_transcript ,level)) {
				matching_pred_transcripts.add(candidate_transcript);
			}
		}
		return(matching_pred_transcripts);
	}
	
	private boolean does_transcript_match(String tr_id, Level level) {
		//System.out.println("Transcript " + transcript_id + "\npredicted transcript " + tr_id);
		if (exons.size() != existing_transcripts.get(tr_id).getExons().size()) { 
			return false;
		}
		List<GTFRecord> tr_exons = existing_transcripts.get(tr_id).getExons();
		for (int e = 0; e < exons.size(); e++) {
			if (e == 0 && level == Level.FLEXIBLE) {
				if ((exons.get(e).strand == Strand.Strand_Positive && exons.get(e).end_pos != tr_exons.get(e).end_pos) ||
						(exons.get(e).strand == Strand.Strand_Negative && exons.get(e).end_pos != tr_exons.get(e).end_pos) ) {
					return false;
				}
			}
			else if (e == exons.size() - 1 && level == Level.FLEXIBLE) {
				if ((exons.get(e).strand == Strand.Strand_Positive && exons.get(e).start_pos != tr_exons.get(e).start_pos) || 
						(exons.get(e).strand == Strand.Strand_Negative && exons.get(e).start_pos != tr_exons.get(e).start_pos)) {
					return false;
				}
			}
			else {
				if (exons.get(e).start_pos != tr_exons.get(e).start_pos || exons.get(e).end_pos != tr_exons.get(e).end_pos) {
					return false;
				}
			}
		}
		return true;
	}
	
	// TODO container = gene-id
	public GTFTranscript(GTFRecord exon) {
		super(exon.sequence, RecordType.TRANSCRIPT.toString(), exon.start_pos, exon.end_pos, exon.strand.toString(), null);
		exons = new ArrayList<GTFRecord>();
		exons.add(exon);
	}
	
	public GTFTranscript(GTFRecord exon, String gene_id) {
		super(exon.sequence, RecordType.TRANSCRIPT.toString(), exon.start_pos, exon.end_pos, exon.strand.toString(), gene_id);
		exons = new ArrayList<GTFRecord>();
		exons.add(exon);
	}
	
	public void addExon(GTFRecord exon) {
		exons.add(exon);
		if (exon.end_pos > this.end_pos) {
			this.end_pos = exon.end_pos;
		}
		if (exon.start_pos < this.start_pos) {
			this.start_pos = exon.start_pos;
		}
	}
	
	public List<GTFRecord> getExons() {
		return this.exons;
	}
	
	public int getLength() {
		int length = 0;
		for (GTFRecord exon : exons) {
			length += exon.getLength();
		}
		return length;
	}
}
