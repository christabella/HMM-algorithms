import sys

'''
Given the current state probability distribution, return the probabilities for the 
different emissions after the next transition.
'''


def hmm1():
    # Parse space-separated stdin values into arrays
    transition_line = sys.stdin.readline().split()
    emission_line = sys.stdin.readline().split()
    pi_line = sys.stdin.readline().split()

    # Reshape into needed matrices
    transition_matrix = reshape(transition_line[2:], int(transition_line[0]), int(transition_line[1]))
    emission_matrix = reshape(emission_line[2:], int(emission_line[0]), int(emission_line[1]))
    pi_matrix = reshape(pi_line[2:], int(pi_line[0]), int(pi_line[1]))

    # Multiply
    pi_A = multiply_matrices(pi_matrix, transition_matrix)
    pi_A_B = multiply_matrices(pi_A, emission_matrix)

    result = str(len(pi_A_B)) + " " + str(len(pi_A_B[0])) + " " # rows and cols
    flattened_pi_A_B = [str(round(elem, 6)) for row in pi_A_B for elem in row]
    result += ' '.join(flattened_pi_A_B)

    return result

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
    elems = [float(x) for x in elems]
    return [elems[i:i+cols] for i in range(0, len(elems), cols)]

def multiply_matrices(A, B):
    """
    Multiply two 2-dimensional matrices.
    """
    return [[sum(a*b for a,b in zip(A_row,B_col)) for B_col in zip(*B)] for A_row in A]

print(hmm1())