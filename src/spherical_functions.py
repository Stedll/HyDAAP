def get_coeffs(signal, filename, bandwidth=200):
    import numpy as np
    import pyshtools as pysh
    grid = pysh.SHGrid.from_array(signal)
    coeffs = grid.expand(lmax_calc=bandwidth-1, nthreads=0)

    array_real = np.zeros(bandwidth*bandwidth)
    array_imag = np.zeros(bandwidth*bandwidth)

    array_real[:int(bandwidth*(bandwidth+1)/2)] = np.split(coeffs.to_array(),2)[0].T[np.triu_indices(bandwidth)].reshape(-1)
    array_imag[:int(bandwidth*(bandwidth+1)/2)] = np.split(coeffs.to_array(),2)[1].T[np.triu_indices(bandwidth)].reshape(-1)
    index = lambda m, l, bandwidth: (m * (bandwidth) - ((m * (m - 1)) / 2) + (l - m)) if m >= 0 else ((((bandwidth - 1) * (bandwidth + 2)) / 2) + 1 + (((bandwidth - 1 + m) * (bandwidth + m) / 2) + (l - abs(m))))
    index_1 = np.asarray([index(m,l,bandwidth) for m in range(1, bandwidth) for l in range(m, bandwidth)], dtype=int)
    index_2 = np.asarray([index(-m,l,bandwidth) for m in range(1, bandwidth) for l in range(m, bandwidth)], dtype=int)
    pow =     np.asarray([1-2*(2%2) for m in range(1, bandwidth) for l in range(m, bandwidth)], dtype=int)

    array_real[index_2] = pow * array_real[index_1];
    array_imag[index_2] = -pow * array_imag[index_1];

    data = np.empty((array_real.size + array_imag.size,), dtype=array_real.dtype)
    data[0::2] = array_real
    data[1::2] = array_imag
    np.savetxt(filename, data, fmt='%.8e', delimiter='\n')
    
    return filename

def rotate_coeffs(signal, euler, body):
    import numpy as np
    import pyshtools as pysh
    grid = pysh.SHGrid.from_array(signal)
    coeffs = grid.expand(nthreads=0)
    coeffs_rotated = coeffs.rotate(euler[0], euler[1], euler[2], degrees=False, convention='y', body=body, dj_matrix=pysh.rotate.djpi2(coeffs.lmax), nthreads=0)
    #coeffs_rotated = pysh.SHCoeffs.from_array(pysh.rotate.SHRotateRealCoef(coeffs.to_array(), euler, pysh.rotate.djpi2(coeffs.lmax)))
    grid = coeffs_rotated.expand(nthreads=0)

    return grid.to_array()

if __name__ == '__main__':
    print('code not implemented for standalone functionalities')
    quit()