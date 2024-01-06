#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

using namespace std;

const size_t MAX_SEQ_LEN = 10000;

// https://github.com/guilhermeagostinelli/levenshtein/blob/master/levenshtein.cpp
// Returns the Levenshtein distance between word1 and word2. 
int levenshteinDist(string &word1, string &word2) {
    int size1 = word1.size();
    int size2 = word2.size();
    int verif[size1 + 1][size2 + 1]; // Verification matrix i.e. 2D array which will store the calculated distance.

    // If one of the words has zero length, the distance is equal to the size of the other word.
    if (size1 == 0)
        return size2;
    if (size2 == 0)
        return size1;

    // Sets the first row and the first column of the verification matrix with the numerical order from 0 to the length of each word.
    for (int i = 0; i <= size1; i++)
        verif[i][0] = i;
    for (int j = 0; j <= size2; j++)
        verif[0][j] = j;

    // Verification step / matrix filling.
    for (int i = 1; i <= size1; i++) {
        for (int j = 1; j <= size2; j++) {
            // Sets the modification cost.
            // 0 means no modification (i.e. equal letters) and 1 means that a modification is needed (i.e. unequal letters).
            int cost = (word2[j - 1] == word1[i - 1]) ? 0 : 1;

            // Sets the current position of the matrix as the minimum value between a (deletion), b (insertion) and c (substitution).
            // a = the upper adjacent value plus 1: verif[i - 1][j] + 1
            // b = the left adjacent value plus 1: verif[i][j - 1] + 1
            // c = the upper left adjacent value plus the modification cost: verif[i - 1][j - 1] + cost
            verif[i][j] = min(
                min(verif[i - 1][j] + 1, verif[i][j - 1] + 1),
                verif[i - 1][j - 1] + cost
            );
        }
    }

    // The last position of the matrix will contain the Levenshtein distance.
    return verif[size1][size2];
}

bool condition(size_t len1, size_t len2, size_t dif) {
    return dif <= 0.1*min(len1, len2);
}

int main() {
    ifstream file("file.fasta");
    vector<string> sequences;
    string line, sequence;
    while (getline(file, line)) {
        if (line[0] == '>') {
            if (!sequence.empty()) {
                sequences.push_back(sequence);
                sequence.clear();
            }
        } else {
            sequence += line;
        }
    }
    if (!sequence.empty()) {
        sequences.push_back(sequence);
    }
    file.close();
    for (auto cur_seq : sequences) {
        cout << cur_seq << '\n';
        size_t cur_seq_len = cur_seq.size();
        size_t max_len = min(cur_seq_len, MAX_SEQ_LEN);
        for (size_t len = max_len; len >= 2; len--) {
            size_t mid = len/2;
            size_t deviation = ceil(0.05*len);
            for (size_t start_pos = 0; start_pos <= cur_seq_len - len; start_pos++) {
                bool decrease = true;
                bool foundAns = false;
                size_t border_point = start_pos + mid;
                for (size_t right_pos = mid - deviation; right_pos < mid + deviation; right_pos++) {
                    size_t len1 = right_pos - start_pos;
                    size_t len2 = len - len1; 
                    string s1 = cur_seq.substr(start_pos, len1);
                    string s2 = cur_seq.substr(right_pos, len2);
                    size_t lev_dist = levenshteinDist(s1, s2);
                    cout << s1 << ' ' << s2 << ' ' << lev_dist << ' ' << deviation << '\n';
                    if (condition(len1, len2, lev_dist)) {
                        foundAns = true;
                        cout << start_pos << ' ' << len1 << ' ' << len2 << '\n';
                        break;
                    }
                }
                if (foundAns)
                    break;
            }
        }
    }
    return 0;
}