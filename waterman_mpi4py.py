from mpi4py import MPI

# Function to compute the alignment score for two characters
def get_alignment_score(a, b):
    # Return a positive or negative score based on the match/mismatch
    match_score = 1
    mismatch_score = -1
    return match_score if a == b else mismatch_score

# Function to compute the Waterman DNA alignment score for a portion of the sequences
def waterman_dna(seq1, seq2, start, end):
    len1 = len(seq1)
    len2 = len(seq2)

    # Create a 2D matrix to store alignment scores
    dp = [[0 for _ in range(len2 + 1)] for _ in range(len1 + 1)]

    # Fill in the matrix using dynamic programming for the given portion of sequences
    for i in range(start + 1, end + 1):
        for j in range(1, len2 + 1):
            match_score = dp[i - 1][j - 1] + get_alignment_score(seq1[i - 1], seq2[j - 1])
            delete_score = dp[i - 1][j]
            insert_score = dp[i][j - 1]

            dp[i][j] = max(0, match_score, delete_score, insert_score)

    # Find the maximum score in the portion of the matrix
    max_score = max(dp[i][j] for i in range(start + 1, end + 1) for j in range(len2 + 1))

    # Reduce the local maximum scores to find the global maximum score within the MPI process
    local_max_score = comm.reduce(max_score, op=MPI.MAX, root=0)

    return local_max_score

# Function to parallelize the Waterman DNA algorithm using MPI
def parallel_waterman_dna(seq1, seq2):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    len1 = len(seq1)
    len2 = len(seq2)
    max_score = 0

    # Distribute workload across nodes using scatter
    local_size = len1 // size
    start = rank * local_size
    end = start + local_size - 1 if rank != size - 1 else len1 - 1

    # Compute the local maximum score for each node
    local_max_score = waterman_dna(seq1, seq2, start, end)

    # Reduce the local maximum scores to find the global maximum score
    max_score = comm.reduce(local_max_score, op=MPI.MAX, root=0)

    return max_score

if __name__ == "__main__":
    sequence1 = "ATCGATCGATCG"
    sequence2 = "CGATCGATCGAT"

    # MPI Initialization
    comm = MPI.COMM_WORLD

    result = parallel_waterman_dna(sequence1, sequence2)

    # Print the result only on the root process (rank 0)
    if comm.rank == 0:
        print("Alignment Score:", result)
