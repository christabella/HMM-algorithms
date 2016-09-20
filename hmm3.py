import sys

"""
Viterbi algorithm
"""


def hmm3():
    # Parse space-separated stdin values into arrays of floats
    transition_line = [float(x) for x in sys.stdin.readline().split()]
    emission_line = [float(x) for x in sys.stdin.readline().split()]
    pi_line = [float(x) for x in sys.stdin.readline().split()]

    emission_sequence = [int(x) for x in sys.stdin.readline().split()[1:]]

    # Reshape into needed matrices
    A = reshape(transition_line[2:], int(transition_line[0]), int(transition_line[1]))
    B = reshape(emission_line[2:], int(emission_line[0]), int(emission_line[1]))
    pi_matrix = reshape(pi_line[2:], int(pi_line[0]), int(pi_line[1]))

    state_count = len(A)


    def viterbi_recursive(delta, delta_idx, emission_seq):
        """
        delta : probabilities 
        delta_idx : states

        all_probabilities : 
        [[0.0, 0.0, 0.0, 0.0], [0.07200000000000002, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]

        """
        # not done... brain hurts...

        if len(emission_seq) == 0:
            # backtrack
            state_sequence = []

            last_state = delta.index(max(delta))
            state_sequence.append(last_state)
            for i in range(len(delta_idx) - 1, 0, -1):
                state_sequence.insert(0, delta_idx[i][last_state])
                last_state = delta_idx[i][last_state]

            print(' '.join([str(x) for x in state_sequence]))
            return

        # implementing max_j∈[1,..N] a_j,i * δt−1(j) * b_i(ot)
        all_probabilities = [[delta[prev_st] * A[prev_st][curr_st] * B[curr_st][emission_seq[0]] for prev_st in range(state_count)] for curr_st in range(state_count)]
        max_probability = [max(probabilities_curr_st) for probabilities_curr_st in all_probabilities]
        # In order to be able to trace back the most likely sequence later on, it is convenient to store the indices
        # of the most likely states at each step.
        delta_idx.append([probabilities_curr_st.index(max_probability[i]) for i, probabilities_curr_st in enumerate(all_probabilities)])
        
        viterbi_recursive(max_probability, delta_idx, emission_seq[1:])

    initial_delta = multiply(pi_matrix[0], get_column(B, emission_sequence[0]))
    initial_delta_idx = [[None]*state_count]

    # Implement Viterbi algorithm
    viterbi_recursive(initial_delta, initial_delta_idx, emission_sequence[1:])

def reshape(elems, rows, cols):
    """
    Reshape a 1-dimensional list into a nested 2-dimensional list.

    Parameters
    ----------
    elems : list
        List of elements to reshape into a 2-D matrix 
    rows : int
        Number of rows in reshaped matrix.
    cols : int
        Number of columns in reshaped matrix.
    """
    return [elems[i:i+cols] for i in range(0, len(elems), cols)]

def multiply(A, B):
    """ 
    Multiply two lists, element-wise, returning a list.
    """
    return [a*b for a, b in zip(A, B)]

def get_column(matrix, col):
    return [row[col] for row in matrix]

hmm3()