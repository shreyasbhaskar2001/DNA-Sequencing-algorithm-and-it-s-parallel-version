import numpy as np
from numba import njit, prange

# Function to compute the alignment score for two characters
@njit
def get_alignment_score(a, b):
    # Return a positive or negative score based on the match/mismatch
    match_score = 1
    mismatch_score = -1
    return match_score if a == b else mismatch_score

# Function to compute the Waterman DNA alignment score for a portion of the sequences
@njit(parallel=True)
def waterman_dna(seq1, seq2, start, end):
    len1 = len(seq1)
    len2 = len(seq2)

    # Create a 2D matrix to store alignment scores
    dp = np.zeros((len1 + 1, len2 + 1), dtype=np.int32)

    # Fill in the matrix using dynamic programming for the given portion of sequences
    for i in prange(start + 1, end + 1):
        for j in range(1, len2 + 1):
            match_score = dp[i - 1, j - 1] + get_alignment_score(seq1[i - 1], seq2[j - 1])
            delete_score = dp[i - 1, j]
            insert_score = dp[i, j - 1]

            dp[i, j] = max(0, match_score, delete_score, insert_score)

    # Find the maximum score in the portion of the matrix
    max_score = np.max(dp[start + 1:end + 1, :])

    return max_score

# Function to parallelize the Waterman DNA algorithm using Numba
def parallel_waterman_dna(seq1, seq2):
    len1 = len(seq1)
    len2 = len(seq2)
    max_score = 0

    # Distribute workload across CPU threads
    chunk_size = len1 // 4  # You can adjust this based on your CPU's core count
    starts = list(range(0, len1, chunk_size))
    ends = [min(start + chunk_size, len1) for start in starts]

    # Compute the local maximum score for each thread
    max_scores = [waterman_dna(seq1, seq2, starts[i], ends[i]) for i in range(len(starts))]

    # Find the global maximum score
    max_score = max(max_scores)

    return max_score

if __name__ == "__main__":
    sequence1 = "ATCGATCGATCG"
    sequence2 = "CGATCGATCGAT"

    result = parallel_waterman_dna(sequence1, sequence2)
    print("Alignment Score:", result)
