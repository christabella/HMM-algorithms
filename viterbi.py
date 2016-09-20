import pprint

def viterbi(obs, states, pi_matrix, transition_matrix, emission_matrix):
    '''
    V is an array of dicts where each dict stores  
        [value] : the previous state and probability of going from the previous state to the current state
        [key]   : the current state

    e.g. 
        [{0: {'prev': None, 'prob': 0.1},
          1: {'prev': None, 'prob': 0.0},
          2: {'prev': None, 'prob': 0.0},
          3: {'prev': None, 'prob': 0.0}},
         {0: {'prev': 0, 'prob': 0.0},
          1: {'prev': 0, 'prob': 0.07200000000000002},
          2: {'prev': 0, 'prob': 0.0},
          3: {'prev': 0, 'prob': 0.0}},
         {0: {'prev': 1, 'prob': 0.0},
          1: {'prev': 0, 'prob': 0.0},
          2: {'prev': 1, 'prob': 0.05184000000000002},
          3: {'prev': 1, 'prob': 0.0}}
        ]
        '''
    V = [{}]
    # Initialize V
    for st in states:
        V[0][st] = {"prob": pi_matrix[st] * emission_matrix[st][obs[0]], "prev": None}
    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
        V.append({})
        # Find maximum probability going from previous state to next
        print("\n")
        for st in states:
            # >> V[t-1][prev_st]["prob"]*transition_matrix[prev_st][st] for prev_st in states
            # In plain writing:
            # Multiply previously calculated probability by probability of going from that state to this 
            max_tr_prob = max(V[t-1][prev_st]["prob"]*transition_matrix[prev_st][st] for prev_st in states)
            print([V[t-1][prev_st]["prob"]*transition_matrix[prev_st][st] for prev_st in states])
            for prev_st in states:
                if V[t-1][prev_st]["prob"]*transition_matrix[prev_st][st] == max_tr_prob:
                    max_prob = max_tr_prob * emission_matrix[st][obs[t]]
                    V[t][st] = {"prob": max_prob, "prev": prev_st}
                    break
    # for line in dptable(V):
    #     print line
    pprint.pprint(V)
    opt = []
    # The highest probability   
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None
    # Get most probable state and its backtrack
    for st, data in V[-1].items():
        if data["prob"] == max_prob:
            opt.append(st)
            previous = st
            break
    # Follow the backtrack till the first observation
    for t in range(len(V) - 2, -1, -1):
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]

    print 'The steps of states are ' + ' '.join([str(x) for x in opt]) + ' with highest probability of %s' % max_prob

def dptable(V):
    # Print a table of steps from dictionary
    yield " ".join(("%12d" % i) for i in range(len(V)))
    for state in V[0]:
        yield "%.7s: " % state + " ".join("%.7s" % ("%f" % v[state]["prob"]) for v in V)

obs = [1, 1, 2, 2]
states = [0, 1, 2, 3]
transition_matrix = [[0.0, 0.8, 0.1, 0.1], [0.1, 0.0, 0.8, 0.1], [0.1, 0.1, 0.0, 0.8], [0.8, 0.1, 0.1, 0.0]]
emission_matrix = [[0.9, 0.1, 0.0, 0.0], [0.0, 0.9, 0.1, 0.0], [0.0, 0.0, 0.9, 0.1], [0.1, 0.0, 0.0, 0.9]]
pi_matrix = [1.0, 0.0, 0.0, 0.0]

viterbi(obs, states, pi_matrix, transition_matrix, emission_matrix)