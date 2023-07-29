#include <iostream>
#include <vector>
#include <string>
#include <omp.h>

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

// Function to compute the Waterman DNA alignment score for two sequences
int watermanDNA(const std::string& seq1, const std::string& seq2) {
    int len1 = seq1.length();
    int len2 = seq2.length();

    // Create a 2D matrix to store alignment scores
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1, 0));

    // Fill in the matrix using dynamic programming
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            int match_score = dp[i - 1][j - 1] + getAlignmentScore(seq1[i - 1], seq2[j - 1]);
            int delete_score = dp[i - 1][j];
            int insert_score = dp[i][j - 1];

            dp[i][j] = std::max(0, std::max(match_score, std::max(delete_score, insert_score)));
        }
    }

    // Find the maximum score in the matrix
    int max_score = 0;
    for (int i = 0; i <= len1; i++) {
        for (int j = 0; j <= len2; j++) {
            if (dp[i][j] > max_score) {
                max_score = dp[i][j];
            }
        }
    }

    return max_score;
}

// Function to parallelize the Waterman DNA algorithm
int parallelWatermanDNA(const std::string& seq1, const std::string& seq2) {
    int len1 = seq1.length();
    int len2 = seq2.length();
    int max_score = 0;

    // Outer loop parallelization (parallelize for different starting points of seq1)
    #pragma omp parallel for shared(max_score)
    for (int i = 0; i < len1; i++) {
        int local_max_score = 0;

        // Inner loop parallelization (parallelize for different starting points of seq2)
        #pragma omp parallel for reduction(max:local_max_score)
        for (int j = 0; j < len2; j++) {
            // Calculate the score for the alignment starting at (i, j)
            int score = watermanDNA(seq1.substr(i), seq2.substr(j));
            
            // Update the local maximum score
            if (score > local_max_score) {
                local_max_score = score;
            }
        }

        // Update the global maximum score
        #pragma omp critical
        {
            if (local_max_score > max_score) {
                max_score = local_max_score;
            }
        }
    }

    return max_score;
}

int main() {
    std::string sequence1 = "ATCGATCGATCG";
    std::string sequence2 = "CGATCGATCGAT";

    int result = parallelWatermanDNA(sequence1, sequence2);
    std::cout << "Alignment Score: " << result << std::endl;

    return 0;
}
