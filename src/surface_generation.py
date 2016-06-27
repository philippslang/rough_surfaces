import numpy as np
import surface as sr


def self_affine():
    '''
    edge_length, roughness_RMS, hurst, N_power_of_two, seed=None, lambda_L_over_lambda_0 =0, lambda_L_over_lambda_1=0, lambda_L_over_lambda_mm=0
    '''
    use_calgo = True
    # defaults
    if seed != None:
        np.random.seed(seed)
    else:
        seed = 0 # tells calgo to use system time (~random seed)
    lambda_L_over_lambda_0 = 1 if lambda_L_over_lambda_0 < 1 else lambda_L_over_lambda_0
    lambda_L_over_lambda_1 = sys.maxint if lambda_L_over_lambda_1 < 1 else lambda_L_over_lambda_1
    mismatched = True if lambda_L_over_lambda_mm > 0 else False
    # lifting
    N = 2**N_power_of_two        
    power = -(hurst+1.)    
    q_L = 2*np.pi/edge_length # q = abs freq
    q_mm = q_L*lambda_L_over_lambda_mm
    f_L, f_d = 1./N, np.sqrt(2*0.5**2) # f = rel fre
    f_0, f_1 = f_L*lambda_L_over_lambda_0, f_L*lambda_L_over_lambda_1        
    if use_calgo:
        A_real, A_imag = np.zeros((N,N)), np.zeros((N,N))
        calgo.self_affine_psd_based_ext(A_real, A_imag, power, seed, f_L, q_L, f_0, f_1, mismatched, q_mm)
        A = A_real + A_imag*1j
    else:
        A = np.zeros((N,N), dtype=complex)
        rand_norm_1, rand_norm_2 = np.random.randn(N/2+1,N/2+1), np.random.randn(N/2+1,N/2+1)
        rand_unif_1, rand_unif_2 = np.random.rand(N/2+1,N/2+1), np.random.rand(N/2+1,N/2+1)
        if mismatched:
            rand_unif_m1, rand_unif_m2 = np.random.rand(N/2+1,N/2+1), np.random.rand(N/2+1,N/2+1)
        for i in range(0, N/2+1):
            for j in range(0, N/2+1):
                phase = 2. * np.pi * rand_unif_1[i,j]                       
                rad = 0.
                f = np.sqrt((float(i)/N)**2 + (float(j)/N)**2)
                if i != 0 or j != 0:       
                    f = f if f > f_0 else f_0 # hi pass for f_0 power
                    rad = rand_norm_1[i,j] * f**power
                if f > f_1: # lo pass for f_1 
                    rad, phase = 0., 0.
                if mismatched: # lo pass for f_mm
                    phase += 2.*np.pi*rand_unif_m1[i,j]*_gamma(f*q_L/f_L, q_mm)
                A[i , j] = rad*np.cos(phase) + rad*np.sin(phase)*1j
                i0 = 0 if i == 0 else N-i
                j0 = 0 if j == 0 else N-j
                A[i0,j0] = rad*np.cos(phase) - rad*np.sin(phase)*1j 
        A[N/2,0] = A[N/2,0].real + 0j
        A[0,N/2] = A[0,N/2].real + 0j
        A[N/2,N/2] = A[N/2,N/2].real + 0j
        for i in range(1, N/2):
            for j in range(1, N/2):
                phase = 2. * np.pi * rand_unif_2[i,j]
                f = np.sqrt((float(i)/N)**2 + (float(j)/N)**2)
                f = f if f > f_0 else f_0 # hi pass for f_0 power            
                rad = rand_norm_2[i,j] * f**power
                if f > f_1:  # lo pass for f_1 
                    rad, phase = 0., 0.
                if mismatched: # lo pass for f_mm
                    phase += 2.*np.pi*rand_unif_m2[i,j]*_gamma(f*q_L/f_L, q_mm)
                A[i,N-j] = rad*np.cos(phase) + rad*np.sin(phase)*1j
                A[N-i,j] = rad*np.cos(phase) - rad*np.sin(phase)*1j
    H = np.real(np.fft.ifft2((A)))
    h_grid = sm.SurfaceMeshSquareRegular(N, edge_length/float(N))
    h_grid.h = H
    h_grid.scale_to_rms(roughness_RMS)
    h_grid.normalize_h()
    return h_grid


if __name__ == '__main__':
    import doctest
    doctest.testmod()