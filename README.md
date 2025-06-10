# Simple Sequence Analyzer

## Project Description

A Python module for basic genomic sequence analysis and visualization. Biological sequences are represented as objects, allowing for intuitive manipulation and analysis. Key functionalities include calculating GC content, performing reverse complement, transcription, and integrating with NCBI's BLAST service for sequence similarity searches.

## Features

* **OOP Design:** Classes for generic sequences, DNA sequences, and genes, demonstrating a clear class hierarchy.
* **Core Bioinformatics Analysis:** Calculate sequence length, GC content, perform transcription (DNA to RNA), and generate reverse complements.
* **BLAST Integration:** Perform online BLASTN similarity searches against NCBI nucleotide databases (`nt`).
* **Basic Visualization:** Plot GC content over sequence windows using `matplotlib` and visualize gene locations within a sequence using ASCII art.

## Coding Features & Best Practices

* **Clear Class Hierarchy:** Demonstrates the use of a base class (`Sequence`) and specialized subclasses (`DNASequence`, `Gene`), illustrating polymorphism and inheritance.
* **Encapsulation:** Each class encapsulates its data (attributes) and behavior (methods), promoting modular and maintainable code.
* **Fundamental Bioinformatics Concepts:** Covers essential operations like sequence representation, length, GC content, reverse complement, and transcription.
* **Simple Visualization:** Shows basic data visualization capabilities using `matplotlib` (for numerical data) and ASCII art (for sequence features).
* **Well-Organized:** The file structure promotes clarity and reusability, forming a proper Python package.
* **Extensibility:** The design allows for easy addition of new sequence types (e.g., `RNASequence`, `ProteinSequence`) or analysis methods in the future without significant refactoring.

## Installation & Setup

To get a local copy up and running, follow these simple steps.

### Prerequisites

* Python 3.8+
* `pip` (Python package installer)
* `git` (for cloning the repository)

### Steps to Install and Run

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/your-username/simple-sequence-analyzer.git](https://github.com/your-username/simple-sequence-analyzer.git) # Replace with your actual repo URL
    cd simple_sequence_analyzer
    ```

2.  **Create and Activate a Virtual Environment:**
    It's highly recommended to use a virtual environment to manage dependencies and avoid conflicts with other Python projects.

    ```bash
    python -m venv venv
    ```
    * **On Windows (Git Bash/MinGW64):**
        ```bash
        source venv/Scripts/activate
        ```
    * **On macOS/Linux:**
        ```bash
        source venv/bin/activate
        ```
    Your terminal prompt should now show `(venv)` indicating the virtual environment is active.

3.  **Install Dependencies:**
    Once your virtual environment is active, install the required packages using `pip`:
    ```bash
    pip install -r requirements.txt
    ```
    If you don't have a `requirements.txt` yet, create one by running `pip freeze > requirements.txt` after installing `biopython` and `matplotlib`:
    ```bash
    pip install biopython matplotlib
    pip freeze > requirements.txt
    ```
    (Note: `pip install --upgrade pip` is generally good practice but doesn't need to be in the README unless it's a specific fix).

## Usage

This module can be imported into your Python scripts or run directly via `examples.py` to see its functionalities in action.

### Running the Example Script

Navigate to the **project root directory** (the `simple_sequence_analyzer` directory containing `examples.py` and the `simple_sequence_analyzer/` package folder) with your virtual environment activated, then run:

```bash
(venv) python -m simple_sequence_analyzer.examples
```

# Usage - Example Code Snippets

Here's how you can use the `simple_sequence_analyzer` module:

```python
from simple_sequence_analyzer.sequences import DNASequence, Gene
from simple_sequence_analyzer.visualization import plot_gc_content_window, visualize_gene_on_sequence

# Create a DNA Sequence object
my_dna = DNASequence("MyGenome", "ATGCGTACGTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGCATGCAGCTAGCTAGCTACGTAGCTAGC")
print(my_dna)
print(f"Length: {my_dna.length()}")
print(f"GC Content: {my_dna.gc_content():.2f}%")
print(f"Transcribed RNA: {my_dna.transcribe()[:30]}...") # Showing snippet
print(f"Reverse Complement: {my_dna.reverse_complement()[:30]}...") # Showing snippet

# Create a Gene within the DNA Sequence
my_gene = Gene("Gene_A", my_dna.sequence, 50, 100)
print(f"\nGene: {my_gene.id}")
print(f"Gene Coding Sequence: {my_gene.get_coding_sequence()}")
print(f"Gene GC Content: {my_gene.gc_content():.2f}%")

# Perform a BLAST search (requires internet connection and a sequence likely to hit)
print("\n--- Performing BLAST Search ---")
# Using a known human 18S rRNA fragment that should yield hits
human_18s_rRNA_fragment = DNASequence(
    "Human_18S_rRNA_Fragment",
    "CGTGGGCCCGCCCGGGCTCGCCCGCCGCAGCACCCGCGCCCGCCCGTCGGC"
)
# Relaxing E-value for demonstration to ensure hits are shown
blast_hits = human_18s_rRNA_fragment.blastn_search(database="nt", evalue_threshold=10, limit_hits=5)

if blast_hits:
    print("Significant BLAST hits found:")
    for i, hit in enumerate(blast_hits):
        print(f"  Hit {i+1}: Accession={hit['accession']}, E-value={hit['e_value']:.2e}, Description={hit['description']}")
else:
    print("No significant BLAST hits found or an error occurred.")

# Visualize gene on sequence (ASCII art)
visualize_gene_on_sequence(my_dna, my_gene)

# Plot GC content window (this will open a matplotlib window)
# Uncomment the line below to view the plot:
# plot_gc_content_window(my_dna, window_size=50)
```

## Project Structure
The project is organized into the following directories and files:

```
imple_sequence_analyzer/
├── __init__.py               # Marks the directory as a Python package.
├── sequences.py              # Defines core biological sequence classes:
│                             #   - `Sequence` (base class)
│                             #   - `DNASequence` (inherits from Sequence, adds DNA-specific methods)
│                             #   - `Gene` (inherits from DNASequence, represents a gene region)
├── visualization.py          # Contains functions for visualizing sequence data, e.g.,
│                             #   - `plot_gc_content_window()`
│                             #   - `visualize_gene_on_sequence()`
└── examples.py               # A script demonstrating how to use the module's functionalities.

├── README.md                 # This file, providing project overview and documentation.
├── .gitignore                # Specifies intentionally untracked files to ignore.
└── requirements.txt          # Lists Python dependencies (e.g., Biopython, Matplotlib).
```