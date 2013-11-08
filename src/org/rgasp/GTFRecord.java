package org.rgasp;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GTFRecord implements Comparable<GTFRecord> {
	protected String sequence;
	public enum RecordType {EXON, CDS, TRANSCRIPT, START_CODON, STOP_CODON, GENE};
	protected RecordType record_type;
	public enum RecordClass {FIRST, TERMINAL, INTERNAL, SINGLE};
	protected RecordClass record_class;
	protected int start_pos;
	protected int end_pos;
	public enum Strand{
		Strand_Positive,
		Strand_Negative,
		Strand_Unknown	
	}
	protected Strand strand;
	protected List<String> container;
	protected float overlap = 0.0f;
	
	protected List<GTFRecord> matchingExons;
	protected List<GTFRecord> overlapingExons;

	public GTFRecord(String sequence, String type, int start, int end, String strand, String container) {
		this.sequence = sequence;
		this.record_type = RecordType.valueOf(type);
		this.start_pos = start;
		this.end_pos = end;
		this.strand = Strand.valueOf(strand);
		this.container = new ArrayList<String>();
		this.container.add(container);
		
		matchingExons = new ArrayList<GTFRecord>();
		overlapingExons = new ArrayList<GTFRecord>();
	}
	
	public int getLength() {
		return end_pos - start_pos + 1;
	}
	
	public String getSequence() {
		return sequence;
	}

	public RecordType getRecord_type() {
		return record_type;
	}

	public int getStart_pos() {
		return start_pos;
	}

	public int getEnd_pos() {
		return end_pos;
	}

	public Strand getStrand() {
		return strand;
	}

	public List<String> getContainer() {
		return container;
	}
	
	public void addContainer(String s) {
		container.add(s);
	}

	public int compareTo(GTFRecord other_record) {
		if (this.sequence.compareTo(other_record.sequence) < 0) return -1;
		else if (this.sequence.compareTo(other_record.sequence) > 0) return 1;
		else if (this.strand.compareTo(other_record.strand) < 0) return -1;
		else if (this.strand.compareTo(other_record.strand) > 0) return 1;
		else if (this.start_pos < other_record.start_pos) return -1;
		else if (this.start_pos > other_record.start_pos) return 1;
		else if (this.end_pos < other_record.end_pos) return -1;
		else if (this.end_pos > other_record.end_pos) return 1;
		return 0;
	}
	
	public boolean equals(Object o) {
		if (o instanceof GTFRecord) {
			GTFRecord other_record = (GTFRecord)o;
			return (this.sequence.equals(other_record.sequence) &&
					this.start_pos == other_record.start_pos &&
					this.end_pos == other_record.end_pos &&
					this.strand == other_record.strand);
		}
		else return false;
	}
}
