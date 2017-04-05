
"Smith-Waterman alignment of two sequence, can be DNA or amino acids"
"Also known as viterbi finding the maximum likelihood path of hidden states of a hidden Markov model (HMM)"

"""X, Y are two sequences, S[l1][l2] is substitution scoring,
    g is the gap opening penalty, and ga is the gap extension penalty"""


def align(X, Y, S, g, ga): # do not change this line


    Vm = [[0 for u in range(len(Y))] for t in range(len(X))]
    Rm = [[0 for u in range(len(Y))] for t in range(len(X))]
    Vx = [[0 for u in range(len(Y) + 1)] for t in range(len(X))]
    Rx = [[0 for u in range(len(Y) + 1)] for t in range(len(X))]
    Vy = [[0 for u in range(len(Y))] for t in range(len(X) + 1)]
    Ry = [[0 for u in range(len(Y))] for t in range(len(X) + 1)]

    def SM_function(t, u):
        score = S[X[t]][Y[u]]
        return (score)

    Vm[0][0] = SM_function(0, 0)
    Vx[0][0] = g
    Vy[0][0] = g
    Vm[0][1] = Vy[0][0] + SM_function(0, 1)
    Vm[1][0] = Vx[0][0] + SM_function(1, 0)
    Vx[0][1] = Vy[0][0] + ga
    Vx[1][0] = Vx[0][0] + ga
    Vy[0][1] = Vy[0][0] + ga
    Vy[1][0] = Vx[0][0] + ga

    for t in range(1, len(X)):
        Vx[t][0] = Vx[t - 1][0] + ga

    for u in range(1, len(Y)):
        Vy[0][u] = Vy[0][u - 1] + ga

    for u in range(1, len(Y) + 1):
        Vx[0][u] = Vy[0][u - 1] + ga

    for t in range(1, len(X) + 1):
        Vy[t][0] = Vx[t - 1][0] + ga

    for t in range(1, len(X)):
        Vm[t][0] = Vx[t - 1][0] + SM_function(t, 0)
    for u in range(1, len(Y)):
        Vm[0][u] = Vy[0][u - 1] + SM_function(0, u)

    for t in range(1, len(X)):
        for u in range(1, len(Y)):
            Vm[t][u] = max(Vm[t - 1][u - 1] + SM_function(t, u), Vx[t - 1][u] + SM_function(t, u), Vy[t][u - 1] + SM_function(t, u))
            Vx[t][u] = max(Vm[t - 1][u - 1] + g, Vx[t - 1][u] + ga, Vy[t][u - 1] + ga)
            Vy[t][u] = max(Vm[t - 1][u - 1] + g, Vx[t - 1][u] + ga, Vy[t][u - 1] + ga)

    for t in range(1, len(X)):
        Vx[t][len(Y)] = max(Vm[t - 1][len(Y) - 1] + g, Vx[t - 1][len(Y)] + ga, Vy[t][len(Y) - 1] + ga)

    for u in range(1, len(Y)):
        Vy[len(X)][u] = max(Vm[len(X) - 1][u - 1] + g, Vx[len(X) - 1][u] + ga, Vy[len(X)][u - 1] + ga)

    Rm[0][0] = 'start'
    Rx[0][0] = 'start'
    Ry[0][0] = 'start'

    for t in range(1, len(X)):
        Rm[t][0] = 'x'
        Rx[t][0] = 'x'
        Ry[t][0] = 'x'

    for u in range(1, len(Y)):
        Rm[0][u] = 'y'
        Rx[0][u] = 'y'
        Ry[0][u] = 'y'

    Rx[0][len(Y)] = 'y'
    Ry[len(X)][0] = 'x'

    for t in range(1, len(X)):
        for u in range(1, len(Y)):
            if (Vm[t][u] == Vy[t][u - 1] + SM_function(t, u)):
                Rm[t][u] = 'y'
            else:
                if (Vm[t][u] == Vx[t - 1][u] + SM_function(t, u)):
                    Rm[t][u] = 'x'
                else:
                    if (Vm[t][u] == Vm[t - 1][u - 1] + SM_function(t, u)):
                        Rm[t][u] = 'm'

            if (Vx[t][u] == Vy[t][u - 1] + ga):
                Rx[t][u] = 'y'
            else:
                if (Vx[t][u] == Vx[t - 1][u] + ga):
                    Rx[t][u] = 'x'
                else:
                    if (Vx[t][u] == Vm[t - 1][u - 1] + g):
                        Rx[t][u] = 'm'

            if (Vy[t][u] == Vy[t][u - 1] + ga):
                Ry[t][u] = 'y'
            else:
                if (Vy[t][u] == Vx[t - 1][u] + ga):
                    Ry[t][u] = 'x'
                else:
                    if (Vy[t][u] == Vm[t - 1][u - 1] + g):
                        Ry[t][u] = 'm'

    for t in range(1, len(X)):
        if (Vx[t][len(Y)] == Vy[t][len(Y) - 1] + ga):
            Rx[t][len(Y)] = 'y'
        else:
            if (Vx[t][len(Y)] == Vx[t - 1][len(Y)] + ga):
                Rx[t][len(Y)] = 'x'
            else:
                if (Vx[t][len(Y)] == Vm[t - 1][len(Y) - 1] + g):
                    Rx[t][len(Y)] = 'm'

    for u in range(1, len(Y)):
        if (Vy[len(X)][u] == Vy[len(X)][u - 1] + ga):
            Ry[len(X)][u] = 'y'
        else:
            if (Vy[len(X)][u] == Vx[len(X) - 1][u] + ga):
                Ry[len(X)][u] = 'x'
            else:
                if (Vy[len(X)][u] == Vm[len(X) - 1][u - 1] + g):
                    Ry[len(X)][u] = 'm'

    Vm_end = Vm[len(X) - 1][len(Y) - 1]
    Vx_end = Vx[len(X) - 1][len(Y)]
    Vy_end = Vy[len(X)][len(Y) - 1]

    if (max(Vm_end, Vx_end, Vy_end) == Vy_end):
        end_state = 'y'
    else:
        if (max(Vm_end, Vx_end, Vy_end) == Vx_end):
            end_state = 'x'
        else:
            if (max(Vm_end, Vx_end, Vy_end) == Vm_end):
                end_state = 'm'

    path = []
    if (end_state == 'm'):
        t = len(X) - 1
        u = len(Y) - 1
        path.insert(0, (end_state, t, u, Vm[t][u], X[t], Y[u]))

    if (end_state == 'x'):
        t = len(X) - 1
        u = len(Y)
        path.insert(0, (end_state, t, u, Vx[t][u], X[t], '-'))

    if (end_state == 'y'):
        t = len(X)
        u = len(Y) - 1
        path.insert(0, (end_state, t, u, Vy[t][u], '-', Y[u]))

    current_state = end_state

    while (t > 0 or u > 0):
        if (current_state == 'm'):
            previous_state = Rm[t][u]
            if (previous_state == 'm'):
                t = t - 1
                u = u - 1
                path.insert(0, (previous_state, t, u, Vm[t][u], X[t], Y[u]))

            if (previous_state == 'x'):
                t = t - 1
                path.insert(0, (previous_state, t, u, Vx[t][u], X[t], '-'))

            if (previous_state == 'y'):
                u = u - 1
                path.insert(0, (previous_state, t, u, Vy[t][u], '-', Y[u]))

        if (current_state == 'x'):
            previous_state = Rx[t][u]
            if (previous_state == 'm'):
                t = t - 1
                u = u - 1
                path.insert(0, (previous_state, t, u, Vm[t][u], X[t], Y[u]))

            if (previous_state == 'x'):
                t = t - 1
                path.insert(0, (previous_state, t, u, Vx[t][u], X[t], '-'))

            if (previous_state == 'y'):
                u = u - 1
                path.insert(0, (previous_state, t, u, Vy[t][u], '-', Y[u]))

        if (current_state == 'y'):
            previous_state = Ry[t][u]
            if (previous_state == 'm'):
                t = t - 1
                u = u - 1
                path.insert(0, (previous_state, t, u, Vm[t][u], X[t], Y[u]))

            if (previous_state == 'x'):
                t = t - 1
                path.insert(0, (previous_state, t, u, Vx[t][u], X[t], '-'))

            if (previous_state == 'y'):
                u = u - 1
                path.insert(0, (previous_state, t, u, Vy[t][u], '-', Y[u]))

        current_state = previous_state
    return path
