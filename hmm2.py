import sys
"""
Calculate the probability of observing a certain emission sequence given a HMM model.
"""

def hmm2():
    # Parse space-separated stdin values into arrays of floats
    transition_line = [float(x) for x in sys.stdin.readline().split()]
    emission_line = [float(x) for x in sys.stdin.readline().split()]
    pi_line = [float(x) for x in sys.stdin.readline().split()]

    emission_sequence = [int(x) for x in sys.stdin.readline().split()[1:]]

    # Reshape into needed matrices
    A = reshape(transition_line[2:], int(transition_line[0]), int(transition_line[1]))
    B = reshape(emission_line[2:], int(emission_line[0]), int(emission_line[1]))
    pi_matrix = reshape(pi_line[2:], int(pi_line[0]), int(pi_line[1]))

    def forward_algo(prev_alpha, emission_sequence):
        if len(emission_sequence) == 0:
            print(round(sum(prev_alpha), 6))
            return sum(prev_alpha)
        state_count = len(A)
        # compute [sum_i=1^N alpha_t(i) a_ij]
        first_term = [sum(multiply(prev_alpha, get_column(A, i))) for i in range(state_count)]
        current_alpha = multiply(first_term, get_column(B, emission_sequence[0]))
        forward_algo(current_alpha, emission_sequence[1:])

    # Initialize alpha-pass algorithm
    initial_alpha = multiply(pi_matrix[0], get_column(B, emission_sequence[0]))
    # Implement alpha-pass algorithm
    forward_algo(initial_alpha, emission_sequence[1:])

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

hmm2()