/* Implement pseudocode in Stamp tutorial.
 * 
 * input
 *  A  : transition matrix
 *  B  : observation matrix
 *  pi : initial state distribution
 *  O  : observation sequence
 *  N  : number of states = A_rows, A_cols, B_rows
 *  M  : number of observations = B_cols
 *  T  : length of observation sequence = O_cols
 *
 * output
 *  est_A : estimated transition matrix
 *  est_B : estimated emission matrix
 * 
 */

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>

#include <math.h> /* log*/

using namespace std;

struct Input {
    vector<vector<double>> A;
    vector<vector<double>> B;
    vector<double> pi;
    vector<double> O;
    int N;
    int M;
    int T;
};

struct Output {
    vector<vector<double>> est_A;
    vector<vector<double>> est_B;    
};

Input parse() {
    Input input;
    int N, M, T, temp; 
    string line;
  

    /********** PARSE TRANSITION MATRIX (A) LINE **********/
    getline(cin, line);
    istringstream iss(line);
    iss >> N; // rows of A
    iss >> temp; // cols of A

    vector<vector<double>> A(N, vector<double>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // extract next whitespace-separated element 
            // from iss stream and put it into elem
            iss >> A[i][j];
        }
    }

    /********** PARSE OBSERVATION MATRIX (B) LINE **********/
    getline(cin, line);
    iss.str(line);
    iss >> temp;
    iss >> M; // cols of B

    vector<vector<double>> B(N, vector<double>(M));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            iss >> B[i][j];
        }
    }
    
    /********** PARSE INITIAL STATE DISTRIBUTION (pi) LINE **********/
    getline(cin, line);
    iss.str(line);
    iss >> temp;
    iss >> temp;

    vector<double> pi(N);
    for (int i = 0; i < N; i++) {
        iss >> pi[i]; 
    }

    /********** PARSE OBSERVATION SEQUENCE (O) LINE **********/
    getline(cin, line);
    iss.str(line);
    iss >> T; // length of O

    vector<double> O(T);
    for (int i = 0; i < T; i++) {
        iss >> O[i]; 
    }

    input.A = A;
    input.B = B;
    input.pi = pi;
    input.O = O;
    input.N = N;
    input.M = M;
    input.T = T;
    return input;
}
vector<vector<double>> multiply_matrices(vector<vector<double>> A, vector<vector<double>> B);
vector<double> multiply(vector<double> A, vector<double> B);
vector<double> get_column(vector<vector<double>> A, int col);
void printMatrix(vector<vector<double>> matrix);
void printVector(vector<double> vect);

Output estimate_model(Input input) {
    Output output;
    vector<vector<double>> A = input.A;
    vector<vector<double>> B = input.B;
    vector<double> pi = input.pi;
    vector<double> O = input.O;
    int N = input.N;
    int M = input.M;
    int T = input.T;

    int iters = 0;
    int maxIters = 100;
    double logProb = -99999999999999;
    double oldLogProb = -999999999999999;

    vector<double> c(T);
    vector<vector<double>> alpha(T, vector<double>(N));
    vector<vector<double>> beta(T, vector<double>(N));
    vector<vector<double>> gamma(T, vector<double>(N));
    vector<vector<vector<double>>> di_gamma(T, vector<vector<double>>(N, vector<double>(N)));


    while (iters < maxIters && logProb > oldLogProb) {
        /********** α-PASS **********/
        // Compute α_0
        c[0] = 0;
        vector<double> b_o_0 = get_column(B, O[0]);
        for (int i = 0; i < N; i++) {
            alpha[0][i] = pi[i] * b_o_0[i]; // α_0
            c[0] += alpha[0][i];
        }

        // Scale α_0
        c[0] = 1/c[0];
        for (int i = 0; i < N; i++) {
            alpha[0][i] *= c[0];
        }

        // Compute α_t
        for (int t = 1; t < T; t++) {
            vector<double> b_o_t = get_column(B, O[t]);
            c[t] = 0;
            for (int i = 0; i < N; i++) {
                alpha[t][i] = 0;
                for (int j = 0; j < N; j++) {
                    alpha[t][i] += alpha[t-1][j] * A[j][i];
                }
                alpha[t][i] *= b_o_t[i];
                c[t] += alpha[t][i];
            }
            // Scale α_t
            c[t] = 1/c[t];
            for (int i = 0; i < N; i++) {
                alpha[t][i] *= c[t];
            }
        }

        /********** β-PASS **********/
        // Let βT−1(i) = 1, scaled by cT−1
        for (int i = 0; i < N; i++) {
            beta[T-1][i] = c[T-1];
        }

        // Computer β_t
        for (int t = T-2; t >= 0; t--) {
            vector<double> b_o_tplus1 = get_column(B, O[t+1]);
            for (int i = 0; i < N; i++) {
                beta[t][i] = 0;
                for (int j = 0; j < N; j++) {
                    beta[t][i] += A[i][j] * b_o_tplus1[j] * beta[t+1][j];
                }
                // Scale β_t
                beta[t][i] *= c[t];
            }
        }

        /********** GAMMA AND DI-GAMMA **********/
        double denom;
        // Compute γ_t(i, j) and γ_t(i)
        for (int t = 0; t < T-1; t++) {
            vector<double> b_o_tplus1 = get_column(B, O[t+1]);
            denom = 0;
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    denom += alpha[t][i] * A[i][j] * b_o_tplus1[j] * beta[t+1][j];
                }
            }
            for (int i = 0; i < N; i++) {
                gamma[t][i] = 0;
                for (int j = 0; j < N; j++) {
                    di_gamma[t][i][j] = (alpha[t][i] * A[i][j] * b_o_tplus1[j] * beta[t+1][j])/denom;
                    gamma[t][i] += di_gamma[t][i][j];
                }
            }
        }

        // Special case for γ_T-1(i)
        denom = 0;
        for (int i = 0; i < N; ++i) {
            denom += alpha[T-1][i];
        }
        for (int i = 0; i < N; ++i) {
            gamma[T-1][i] = alpha[T-1][i] / denom;
        }

        /********** Re-estimate A, B and pi⇡ **********/

        // Re-estimate pi
        for (int i = 0; i < N; ++i){
            pi[i] = gamma[0][i];
        }

        double numer;
        // Re-estimate A
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                numer = 0;
                denom = 0;
                for (int t = 0; t < T-1; t++) {
                    numer += di_gamma[t][i][j];
                    denom += gamma[t][i];
                }
                A[i][j] = numer / denom;
            }
        }

        // Re-estimate B
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                numer = 0;
                denom = 0;
                for (int t = 0; t < T; t++) {
                    if (O[t] == j) {
                        numer += gamma[t][i];
                    }
                    denom += gamma[t][i];
                }
                B[i][j] = numer / denom;
            }
        }

        // Compute log[P(O|lambda)]
        oldLogProb = logProb;
        logProb = 0;
        for (int i = 0; i < T; ++i) {
            logProb += log(c[i]);
        }
        logProb = -logProb;

        // To iterate or not to iterate, that is the question...
        iters += 1;
    } // end while

    output.est_A = A;
    output.est_B = B;
    return output;
}

vector<vector<double>> multiply_matrices(vector<vector<double>> A, vector<vector<double>> B) {
    /*
     * ikj-algorithm for matrix multiplication, where
     *  j : rows in A
     *  j : cols in B
     *  k : cols in A = rows in B
     */
    
    int i_size = A.size(); // rows of A and C
    int k_size = B.size(); // rows of B
    int j_size = B[0].size(); // cols of B and C

    vector<vector<double>> C(i_size, vector<double>(j_size));

    for (int i = 0; i < i_size; i++) {
        for (int k = 0; k < k_size; k++) {
            for (int j = 0; j < j_size; j++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

vector<double> multiply(vector<double> A, vector<double> B) {
    /*
     * Multiply two lists, element-wise, returning a list.
     */
    int length = A.size();

    vector<double> C(length);

    for (int i = 0; i < length; i++) {
        C[i] = A[i] * B[i];
    }
    return C;
}

vector<double> get_column(vector<vector<double>> A, int col) {
    int length = A.size();
    vector<double> B(length);
    for (int i = 0; i < length; i++) {
        B[i] = A[i][col];
    }
    return B;
}

void printMatrix(vector<vector<double>> matrix) {
    vector<vector<double>>::iterator it;
    vector<double>::iterator inner;
    for (it=matrix.begin(); it != matrix.end(); it++) {
        for (inner = it->begin(); inner != it->end(); inner++) {
            cout << *inner;
            if(inner+1 != it->end()) {
                cout << "\t";
            }
        }
        cout << endl;
    }
}

void printVector(vector<double> vect) {
    vector<double>::iterator it;
    for (it=vect.begin(); it != vect.end(); it++) {
        cout << *it;
        cout << " ";
    }
}

int main (int argc, char* argv[]) {
    Input input = parse();

    Output output = estimate_model(input);

    cout << input.N;
    cout << ' ';
    cout << input.N;
    for (int i = 0; i < input.N; i++) {
        for (int j = 0; j < input.N; j++) {
            cout << ' ';
            cout << output.est_A[i][j];
        }
    }
    cout << endl;
    cout << input.N;
    cout << ' ';
    cout << input.M;
    for (int i = 0; i < input.N; i++) {
        for (int j = 0; j < input.M; j++) {
            cout << ' ';
            cout << output.est_B[i][j];
        }
    }
    cout << endl;

    // cout << "--------------- Estimated B --------------- \n";
    // printMatrix(output.est_B);
    // cout << endl;

    // cout << "--------------- Estimated A --------------- \n";
    // printMatrix(output.est_A);
    // cout << endl;

    // cout << "--------------- Estimated pi --------------- \n";
    // printVector(input.pi);
    // cout << endl;
    return 0;
}