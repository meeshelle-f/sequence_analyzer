#sequences.py 

from Bio.Blast import NCBIWWW, NCBIXML # For online BLAST queries, parsing XML output from NCBI BLAST
from Bio.Seq import Seq # Biopython's sequence object, good for passing to BLAST
from Bio.SeqRecord import SeqRecord # Useful for BLAST input if you want to include an ID

class Sequence:
    def __init__(self, id: str, sequence: str) -> None:
        """
        Initializes a generic sequence
        
        Args:
            id(str): A unique id for the sequence.
            sequence(str): The actual nucleotdie or amino acid sequence.
        """
        if not isinstance(id, str) or not id:
            raise ValueError("ID must be a non-empty string.")
        if not isinstance(sequence, str) or not sequence:
            raise ValueError("Sequence must be a non-empty string.")
        self.id = id
        self.sequence = sequence.upper()

    def __str__(self):
        return f"ID: {self.id}, Sequence: {self.sequence[:30]}..."  # Truncate for display

    def length(self):
        return len(self.sequence)

    def reverse_complement(self):
        """
        Returns the reverse complement of the sequence.
        This method must be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement reverse_complement.")

class DNASequence(Sequence):
    def __init__(self, id: str, sequence: str) -> None:
        """
        Initializes a DNA sequence.

        Args:
            id (str): A unique identifier for the DNA sequence.
            sequence (str): The DNA sequence string (containing A, T, G, C).
        """     
        super().__init__(id, sequence)
        if not all(nuc in "ATGC" for nuc in self.sequence):
            raise ValueError("DNA sequence can only contain A, T, G, C.")
        self.complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def gc_content(self):
        """
        Calculates and returns the GC content percentage of the DNA sequence.
        """
        g_count = self.sequence.count('G')
        c_count = self.sequence.count('C')
        if not self.length():
            return 0.0
        return ((g_count + c_count) / self.length()) * 100

    def transcribe(self):
        """
        Returns the RNA sequence by replacing 'T' with 'U'.
        """
        return self.sequence.replace('T', 'U')

    def reverse_complement(self):
        """
        Implements the DNA-specific reverse complement.
        """
        complement = "".join([self.complement_map[nuc] for nuc in self.sequence])
        return complement[::-1]

    def blastn_search(self, database: str = "nt", evalue_threshold: float = 0.01, limit_hits: int = 5):
        
        """
        Performs a BLASTN search of this DNA sequence against an NCBI database.

        Args:
            database (str): The NCBI database to search against (e.g., "nt", "nr").
            evalue_threshold (float): The E-value threshold for reported hits.
            limit_hits (int): The maximum number of top hits to return.

        Returns:
            list: A list of dictionaries, where each dictionary represents a BLAST hit
                  with keys like 'accession', 'description', 'e_value', 'bit_score', 'alignment_length'.
                  Returns an empty list if no hits are found or an error occurs.
        """

        print(f"Running BLASTN search for {self.id} against {database} database...")
        try:
            seq_record = SeqRecord(Seq(self.sequence), id=self.id, name=self.id, description=f"Query sequence: {self.id}")

            with NCBIWWW.qblast(
                program="blastn",
                database=database,
                sequence=seq_record.format("fasta"),
                entrez_query=None,
                expect=evalue_threshold,
                hitlist_size=limit_hits,
                format_type="XML"
            ) as result_handle:
                # --- START TEMPORARY DIAGNOSTIC CODE ---
                # Read the entire content of the result_handle
                raw_xml_content = result_handle.read()

                # Print the first 2000 characters (or more if you want)
                print("\n--- RAW BLAST XML CONTENT START ---")
                print(raw_xml_content[:2000]) # Or print(raw_xml_content) for full content
                print("--- RAW BLAST XML CONTENT END ---\n")

                # IMPORTANT: Rewind the handle to the beginning so NCBIXML.parse can read it
                result_handle.seek(0)
                # --- END TEMPORARY DIAGNOSTIC CODE ---

                blast_records = NCBIXML.parse(result_handle)

                hits = []
                for blast_record in blast_records:
                    # If you get here, it means XML parsing *started*.
                    # Check if there are any alignments at all.
                    if not blast_record.alignments:
                        print(f"DEBUG: No alignments found in a blast_record (query: {blast_record.query}).")
                        continue # Move to next blast_record if no alignments

                    for alignment in blast_record.alignments:
                        # Check if there are any HSPs (High-scoring Segment Pairs)
                        if not alignment.hsps:
                            print(f"DEBUG: No HSPs found in alignment (accession: {alignment.accession}).")
                            continue # Move to next alignment if no HSPs

                        for hsp in alignment.hsps:
                            hit_info = {
                                "accession": alignment.accession,
                                "description": alignment.title,
                                "e_value": hsp.expect,
                                "bit_score": hsp.bits,
                                "alignment_length": hsp.align_length,
                            }
                            hits.append(hit_info)
                            if len(hits) >= limit_hits:
                                break
                    if len(hits) >= limit_hits:
                        break
                return hits

        except Exception as e:
            print(f"An error occurred during BLAST search: {e}")
            return []
        
class Gene(DNASequence):
    def __init__(self, id: str, sequence: str, start: int, end: int) -> None:
        """
        Initializes a Gene, representing a specific region of a DNA sequence.

        Args:
            id (str): A unique identifier for the gene.
            sequence (str): The full DNA sequence on which the gene resides.
            start (int): The 0-based start position of the gene (inclusive).
            end (int): The 0-based end position of the gene (exclusive).
        """
        super().__init__(id, sequence)
        if not (0 <= start < end <= self.length()):
            raise ValueError(f"Invalid start/end positions for gene. Sequence length: {self.length()}")
        self.start = start
        self.end = end

    def get_coding_sequence(self) -> str:
        """
        Returns the subsequence corresponding to the gene's coding region.
        """
        return self.sequence[self.start:self.end]
    
   
