#include <iostream>
#include <vector>
#include <string>
#include <mpi.h>

// Function to compute the alignment score for two characters
int getAlignmentScore(char a, char b) {
    const int match_score = 1;
    const int mismatch_score = -1;

    if (a == b) {
        return match_score;
    } else {
        return mismatch_score;
    }
}

// Function to compute the Waterman DNA alignment score for a portion of the sequences
int watermanDNA(const std::string& seq1, const std::string& seq2, int start, int end) {
    int len1 = seq1.length();
    int len2 = seq2.length();

    // Create a 2D matrix to store alignment scores
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1, 0));

    // Fill in the matrix using dynamic programming for the given portion of sequences
    for (int i = start + 1; i <= end; i++) {
        for (int j = 1; j <= len2; j++) {
            int match_score = dp[i - 1][j - 1] + getAlignmentScore(seq1[i - 1], seq2[j - 1]);
            int delete_score = dp[i - 1][j];
            int insert_score = dp[i][j - 1];

            dp[i][j] = std::max(0, std::max(match_score, std::max(delete_score, insert_score)));
        }
    }

    // Find the maximum score in the portion of the matrix
    int max_score = 0;
    for (int i = start + 1; i <= end; i++) {
        for (int j = 0; j <= len2; j++) {
            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
            }
        }
    }

    // Reduce the local maximum scores to find the global maximum score within the MPI process
    int local_max_score = max_score;
    MPI_Allreduce(&local_max_score, &max_score, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    return max_score;
}

// Function to parallelize the Waterman DNA algorithm using MPI
int parallelWatermanDNA(const std::string& seq1, const std::string& seq2) {
    int len1 = seq1.length();
    int len2 = seq2.length();
    int max_score = 0;

    // MPI Initialization
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Distribute workload across nodes using MPI_Scatter
    int local_size = len1 / size;
    int start = rank * local_size;
    int end = (rank == size - 1) ? len1 : start + local_size;

    // Compute the local maximum score for each node
    int local_max_score = watermanDNA(seq1, seq2, start, end);

    // Reduce the local maximum scores to find the global maximum score
    MPI_Reduce(&local_max_score, &max_score, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    // Finalize MPI
    MPI_Finalize();

    return max_score;
}

int main(int argc, char** argv) {
    std::string sequence1 = "ATCGATCGATCG";
    std::string sequence2 = "CGATCGATCGAT";

    // MPI Initialization
    MPI_Init(&argc, &argv);

    int result = parallelWatermanDNA(sequence1, sequence2);

    // Print the result only on the root process (rank 0)
    if (MPI::COMM_WORLD.Get_rank() == 0) {
        std::cout << "Alignment Score: " << result << std::endl;
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
