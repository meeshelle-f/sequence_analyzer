import matplotlib.pyplot as plt
from .sequences import DNASequence, Gene 

def plot_gc_content_window(dna_sequence: DNASequence, window_size: int = 100):
    if not isinstance(dna_sequence, DNASequence):
        raise TypeError("Input must be a DNASequence object.")
    if not isinstance(window_size, int) or window_size <= 0:
        raise ValueError("Window size must be a positive integer.")

    gc_values = []
    for i in range(0, dna_sequence.length() - window_size + 1, window_size // 2): # Sliding window
        window_seq = DNASequence(f"window_{i}", dna_sequence.sequence[i : i + window_size])
        gc_values.append(window_seq.gc_content())

    plt.figure(figsize=(12, 6))
    plt.plot([i for i in range(len(gc_values))], gc_values)
    plt.title(f"GC Content in Sliding Windows for {dna_sequence.id}")
    plt.xlabel(f"Window Index (window size: {window_size}, step: {window_size//2})")
    plt.ylabel("GC Content (%)")
    plt.grid(True)
    plt.show()

def visualize_gene_on_sequence(dna_sequence: DNASequence, gene: Gene): 
    """
    Prints an ASCII representation of a DNA sequence with a gene highlighted.

    Args:
        dna_sequence (DNASequence): The DNASequence object to visualize.
        gene (Union[Gene, str]): A Gene object or a string representing the gene sequence.
    """
        
    if not isinstance(dna_sequence, DNASequence):
        raise TypeError("First argument must be a DNASequence object.")

    gene_sequence: str = ""
    gene_start: int = -1
    gene_end: int = -1

    if isinstance(gene, Gene):
        gene_sequence = gene.get_coding_sequence()
        gene_start = gene.start
        gene_end = gene.end
        gene_id = gene.id

    elif isinstance(gene, str): 
        # This branch handles cases where a raw string is passed,
        # but it will only work if the string is found in the DNA sequence.
        print("Warning: A raw string was provided for 'gene'. Attempting to find its position.")
        gene_sequence = gene
        gene_start = dna_sequence.sequence.find(gene_sequence)
        if gene_start == -1:
            print(f"Error: Gene sequence '{gene_sequence}' not found in the DNA sequence.")
            return # Exit if gene string not found
        gene_end = gene_start + len(gene_sequence)
        gene_id = f"'{gene_sequence}' (raw string)" # Create a temporary ID for display
    else:
        raise TypeError("The 'gene' argument must be a Gene object or a string.")

    print(f"\nSequence: {dna_sequence.id}")
    print(dna_sequence.sequence)

    # Highlight the gene region only if it was found
    if gene_start != -1 and gene_end != -1:
        highlight = ['-'] * dna_sequence.length()
        # Ensure highlight indices are within bounds
        for i in range(max(0, gene_start), min(dna_sequence.length(), gene_end)):
            highlight[i] = '^'
        print("".join(highlight))
        print(f"Gene '{gene_id}' from position {gene_start} to {gene_end}")
    else:
        print("No gene region to highlight.")