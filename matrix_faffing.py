from futils import*

def mat_to_a_power(mat, no_times):
    assert no_times >= 0
    if no_times == 0:
        return np.identity(mat.shape[0], dtype='complex128')
    elif no_times == 1:
        return mat
    else:
        result = mat
        for i in range(no_times-1):
            result = np.dot(result, mat)
        return result


def f_expm(mat, accuracy):
    total = np.zeros(mat.shape, dtype='complex128')
    for i in range(accuracy):
        prefactor = 1 / np.math.factorial(i)
        term = mat_to_a_power(mat, i)
        total += prefactor * term
        print(str(i) + 'th term in expansion is:')
        print(term)
    return total







flipper = Gate('Pauli_x', position=0, n_qubits=2, lead_time=0, out_time=0, duration=0.01)
reg = QRegister(2)

print(flipper.matrix)
print(flipper.gatestring)

reg.output()
flipper.apply(reg)
reg.output()