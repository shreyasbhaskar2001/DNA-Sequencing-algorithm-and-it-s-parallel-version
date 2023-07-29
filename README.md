# DNA-Sequencing-algorithm-and-it-s-parallel-version

Waterman DNA Sequencing Algorithm

The Waterman DNA sequencing algorithm, developed by Michael Waterman in 1976, is a dynamic programming-based approach used to compare and align DNA or protein sequences. It is primarily employed in bioinformatics and computational biology to identify regions of similarity between two sequences, which can provide valuable insights into the evolutionary relationships, functional similarities, and structural properties of biological molecules.

Key Steps of the Waterman DNA Sequencing Algorithm:

Scoring Scheme: The algorithm begins by defining a scoring scheme that assigns scores to matches and mismatches between individual characters (nucleotides or amino acids) in the sequences being compared. The scoring scheme typically assigns positive scores for matches and negative scores for mismatches, with the idea that matches contribute positively to the alignment score while mismatches have a penalizing effect.

Dynamic Programming Matrix: The Waterman algorithm employs a dynamic programming matrix to compute the alignment scores. The matrix is a two-dimensional table with rows and columns corresponding to the length of the two sequences being aligned. The matrix is filled in a systematic manner, considering previously calculated scores to determine the best alignment score for each position.

Filling the Matrix: Starting from the top-left corner of the matrix, the algorithm iteratively fills each cell based on three possible choices: match, deletion, or insertion. The score in each cell is the maximum of the scores obtained from the previous diagonal cell (match/mismatch), the cell above (deletion), and the cell to the left (insertion).

Traceback: After the matrix is fully computed, a traceback step is performed to identify the optimal alignment path that yields the highest alignment score. The traceback starts from the cell with the maximum score and follows the path that led to that score, effectively determining the aligned regions of the sequences.

Alignment Generation: The traceback step generates the aligned sequences by inserting gaps in the appropriate positions to align the characters optimally. The aligned sequences are constructed by walking along the traceback path.

Significance of the Waterman DNA Sequencing Algorithm:

Sequence Comparison: The Waterman algorithm is instrumental in comparing DNA or protein sequences, revealing regions of similarity and dissimilarity, which can help identify conserved functional regions or evolutionary relationships between biological molecules.

Functional Annotation: By identifying conserved regions between known functional sequences and unknown sequences, the algorithm aids in functional annotation, helping to infer the likely function of new DNA or protein sequences.

Structural Biology: In protein studies, sequence alignment is crucial for predicting protein structure and function based on homologous proteins with known structures.

Evolutionary Studies: By aligning homologous sequences from different species, the Waterman algorithm helps in studying evolutionary relationships and constructing phylogenetic trees.

Genome Assembly: The algorithm is a fundamental building block in genome assembly, where it is used to align millions of short DNA reads to reconstruct complete genomes.

Drug Discovery: Sequence alignments are critical in the discovery of drug targets and understanding drug resistance in pathogens.
