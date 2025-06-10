from simple_sequence_analyzer.sequences import DNASequence, Gene
from simple_sequence_analyzer.visualization import plot_gc_content_window, visualize_gene_on_sequence

if __name__ == "__main__":
    # Create a DNA Sequence
    dna_seq = DNASequence("MyGenomicSequence", "ATGCGTACGTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGC")
    print(dna_seq)
    print(f"Length: {dna_seq.length()}")
    print(f"GC Content: {dna_seq.gc_content():.2f}%")
    print(f"Transcribed RNA: {dna_seq.transcribe()[:30]}...")
    print(f"Reverse Complement: {dna_seq.reverse_complement()[:30]}...")

    # Create a Gene within the DNA Sequence
    my_gene = Gene("Gene_A", dna_seq.sequence, 50, 100)
    print(f"\nGene: {my_gene.id}")
    print(f"Gene Coding Sequence: {my_gene.get_coding_sequence()}")
    print(f"Gene GC Content: {my_gene.gc_content():.2f}%")

    # Perform a BLAST search for the DNA sequence
    print("\n--- Performing BLAST Search ---")
    # Use a shorter, more specific sequence for a quicker demo BLAST
    demo_seq = DNASequence("DemoSequence", "ATGCGTAGATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGC")

    test_blast_seq = DNASequence(
        "Human_18S_rRNA_Fragment",
        "CGTGGGCCCGCCCGGGCTCGCCCGCCGCAGCACCCGCGCCCGCCCGTCGGC"
    )

    blast_hits = test_blast_seq.blastn_search(database="nt", evalue_threshold=10, limit_hits=5)


    if blast_hits:
        print("Significant BLAST hits found:")
        for i, hit in enumerate(blast_hits):
            print(f"  Hit {i+1}: Accession={hit['accession']}, E-value={hit['e_value']:.2e}, Description={hit['description']}")
    else:
        print("No significant BLAST hits found.")

    # Visualize GC content
    plot_gc_content_window(dna_seq, window_size=50)

    # Visualize gene on sequence
    visualize_gene_on_sequence(dna_seq, my_gene)

    # Example of a short sequence for gene visualization
    short_dna_seq = DNASequence("ShortSeq", "ATGCATGCATGCATGC")
    short_gene = Gene("ShortGene", short_dna_seq.sequence, 2, 8)
    visualize_gene_on_sequence(short_dna_seq, short_gene)





