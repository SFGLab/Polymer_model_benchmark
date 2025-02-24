import numpy as np
import hicstraw
from tqdm import tqdm
from scipy.spatial import distance
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr, kendalltau
from scipy.ndimage import uniform_filter, gaussian_filter
from sklearn.cross_decomposition import CCA
from scipy.spatial import procrustes
from scipy.stats import wasserstein_distance
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import pairwise_distances
import scipy.spatial.distance as ssd
import scipy.stats as stats
import cv2
from scipy.linalg import sqrtm
from skimage.metrics import structural_similarity as ssim
from scipy.linalg import sqrtm
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy.spatial.distance import pdist, squareform

######### Preprocessing Functions ###########

def get_coordinates_cif(file:str):
    '''
    It returns the corrdinate matrix V (N,3) of a .pdb file.
    The main problem of this function is that coordiantes are not always in 
    the same column position of a .pdb file. Do changes appropriatelly,
    in case that the data aren't stored correctly. 
    
    Input:
    file (str): the path of the .cif file.
    
    Output:
    V (np.array): the matrix of coordinates
    '''
    V = list()
    
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                columns = line.split()
                x = eval(columns[10])
                y = eval(columns[11])
                z = eval(columns[12])
                V.append([x, y, z])
    
    return np.array(V)

def get_heat(V,viz=False):
    '''
    Inputs structure V(Nx3) and outputs the inverse distance heatmap.
    '''
    mat = distance.cdist(V, V, 'euclidean') # this is the way \--/
    mat = 1/(mat+0.1)**(1/3)

    if viz:
        figure(figsize=(15, 7))
        plt.imshow(mat,cmap="Reds",vmax=0.01)
        plt.show()
    return mat

def generate_random_walks(n_walks=1000, n_steps=1000):
    """
    Generate `n_walks` random walk structures in 3D.
    
    Parameters:
    - n_walks (int): Number of random walks to generate.
    - n_steps (int): Number of steps in each random walk.

    Returns:
    - list of numpy arrays: Each array represents a random walk in 3D space.
    """
    walks = []
    heat = 0
    
    for _ in tqdm(range(n_walks)):
        steps = np.random.choice([-1, 1], size=(n_steps, 3))  # Random steps in x, y, z
        positions = np.cumsum(steps, axis=0)  # Compute cumulative positions
        heat += get_heat(positions)
        walks.append(positions)
    
    return walks, heat/n_walks

def remove_diagonal(matrix):
    """
    Removes the diagonal of a square matrix by setting it to zero.

    Parameters:
    - matrix (numpy.ndarray): A square matrix.

    Returns:
    - numpy.ndarray: The matrix with its diagonal elements set to zero.
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Input must be a square matrix.")
    
    matrix_no_diag = matrix.copy()
    np.fill_diagonal(matrix_no_diag, 0)
    
    return matrix_no_diag

def remove_k_diagonals(matrix, k):
    '''
    This function totally removes the diagonal of the heatmap with brute force.
    '''
    N = len(matrix)
    reduced_matrix = np.zeros((N,N))
    reduced_matrix += np.tril(matrix,k=-k)
    reduced_matrix += np.triu(matrix,k=k)
    
    new_matrix = np.zeros((N-k+1,N-k))
    for i in range(k,N):
        for j in range(i-k+1):
            new_matrix[i-k+1,j] = reduced_matrix[i,j]
            new_matrix[j,i-k] = reduced_matrix[j,i]
    new_matrix = new_matrix[:N-k,:N-k]
    return np.array(new_matrix)

def avg_heat_multimm(multimm_ensemble_path,N_ensemble):
    """
    Computes the average distance heatmaps from MultiMM ensemble.
    """
    H = 0
    for i in tqdm(range(N_ensemble)):
        V = get_coordinates_cif(multimm_ensemble_path+f'/sample_{i+1}/MultiMM_minimized.cif')
        H += get_heat(V)
    H = H/N_ensemble
    return H

def avg_heat_loopsage(loopsage_ensemble_path,N_ensemble,mode='EM'):
    """
    Computes the average distance heatmaps from MultiMM ensemble.
    """
    H = 0
    for i in tqdm(range(0,N_ensemble)):
        V = get_coordinates_cif(loopsage_ensemble_path+'/ensemble/'+mode+f'LE_{i+1}.cif')
        H += get_heat(V)
    H = H/N_ensemble
    return H
    

########### Estimate how the signal of the heatmap decreases from the diagonal and derive the power of power law ############

def average_signal_distance(matrix):
    """
    Calculate the average signal as a function of the distance from the diagonal for a symmetric matrix.
    
    Parameters:
        matrix (numpy.ndarray): Symmetric 2D array.
        
    Returns:
        numpy.ndarray: Array containing the average signal for each distance from the diagonal.
    """
    if not isinstance(matrix, np.ndarray) or matrix.ndim != 2:
        raise ValueError("Input must be a 2D numpy array.")
    
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square.")
    
    n = matrix.shape[0]
    distances = np.arange(n)
    signal = np.zeros(n, dtype=float)
    counts = np.zeros(n, dtype=int)
    
    for i in range(n):
        for j in range(i, n):  # Exploit symmetry, only compute for upper triangle
            d = abs(j - i)
            signal[d] += matrix[i, j]
            counts[d] += 1

    # Avoid division by zero
    counts[counts == 0] = 1
    signal /= counts
    return signal

def fit_power_law(distances, p_s):
    """
    Fit a power law p(s) = s^a to the array p_s and estimate the exponent a with error.
    
    Parameters:
        distances (numpy.ndarray): Array of distances from the diagonal (1D).
        p_s (numpy.ndarray): Array of average signals corresponding to distances (1D).
        
    Returns:
        tuple: Estimated exponent `a`, and its standard error.
    """
    if not isinstance(distances, np.ndarray) or not isinstance(p_s, np.ndarray):
        raise ValueError("Both distances and p_s must be numpy arrays.")
    if len(distances) != len(p_s):
        raise ValueError("Distances and p_s must have the same length.")
    if np.any(distances <= 0):
        raise ValueError("Distances must be positive to fit a power law.")
    
    # Power-law function
    def power_law(s, a, c):
        return c * s**a

    # Perform curve fitting
    try:
        popt, pcov = curve_fit(power_law, distances, p_s, p0=[-1.0, 1.0], bounds=([-np.inf, 0], [np.inf, np.inf]))
        a, c = popt
        a_error = np.sqrt(np.diag(pcov))[0]  # Standard error for `a`
    except Exception as e:
        raise RuntimeError(f"Curve fitting failed: {e}")
    print(f'power is {a:.4f} with error +/-{a_error:.4f}')
    return a, a_error

############## Estimation of Regular Correlations between Heatmaps ################

def smooth_matrix(matrix, sigma=1.0):
    """
    Apply Gaussian smoothing to the input matrix.
    
    Parameters:
    - matrix: 2D numpy array to be smoothed.
    - sigma: Standard deviation for Gaussian kernel (controls smoothing extent).
    
    Returns:
    - Smoothed matrix.
    """
    return gaussian_filter(matrix, sigma=sigma)

def normalize_matrix(matrix, method="zscore"):
    """
    Normalize the input matrix using the specified method.
    Methods: 'zscore', 'minmax', 'log'
    """
    if method == "zscore":
        mean = np.mean(matrix)
        std = np.std(matrix)
        return (matrix - mean) / std
    elif method == "minmax":
        min_val = np.min(matrix)
        max_val = np.max(matrix)
        return (matrix - min_val) / (max_val - min_val)
    elif method == "log":
        return np.log1p(matrix)
    else:
        raise ValueError(f"Unknown normalization method: {method}")

def resize_matrix(larger, target_shape):
    """Resize the larger matrix to the target shape using averaging."""
    large_rows, large_cols = larger.shape
    target_rows, target_cols = target_shape

    row_scale = large_rows / target_rows
    col_scale = large_cols / target_cols

    # Create an output matrix
    resized = np.zeros((target_rows, target_cols))

    for i in range(target_rows):
        for j in range(target_cols):
            # Determine the bounds of the block in the original matrix
            row_start = int(i * row_scale)
            row_end = int((i + 1) * row_scale)
            col_start = int(j * col_scale)
            col_end = int((j + 1) * col_scale)

            # Average over the block
            resized[i, j] = np.mean(larger[row_start:row_end, col_start:col_end])

    return resized

def compute_metrics(matrix1, matrix2):
    """Compute various metrics between two matrices."""
    # Flatten the matrices
    m1_flat = matrix1.flatten()
    m2_flat = matrix2.flatten()

    # Pearson correlation
    pearson_corr, _ = pearsonr(m1_flat, m2_flat)

    # Spearman correlation
    spearman_corr, _ = spearmanr(m1_flat, m2_flat)

    # Kendall-Tau correlation
    kendall_corr, _ = kendalltau(m1_flat, m2_flat)

    # Mean Squared Error
    mse = np.mean((m1_flat - m2_flat) ** 2)

    return {
        "Pearson Correlation": pearson_corr,
        "Spearman Correlation": spearman_corr,
        "Kendall-Tau Correlation": kendall_corr,
        "Mean Squared Error": mse,
    }

def compare_heatmaps(hic_exp, hic_sim, sigma=5.0, normalization_method="zscore"):
    """Main function to compare Hi-C experimental and simulated heatmaps."""
    # Apply Gaussian smoothing to the experimental heatmap
    smoothed_exp = smooth_matrix(hic_exp, sigma=sigma)

    # Apply averaging to the bigger matrix to match th dimension of smaller
    if hic_exp.shape == hic_sim.shape:
        resized_exp, resized_sim = smoothed_exp, hic_sim
    elif hic_exp.size > hic_sim.size:
        resized_exp = resize_matrix(smoothed_exp, hic_sim.shape)
        resized_sim = hic_sim
    else:
        resized_sim = resize_matrix(hic_sim, smoothed_exp.shape)
        resized_exp = hic_exp

     # Normalize the matrices
    norm_exp = normalize_matrix(resized_exp, method=normalization_method)
    norm_sim = normalize_matrix(resized_sim, method=normalization_method)

    # Estimate correlation
    metrics = compute_metrics(norm_exp, norm_sim)
    return resized_exp, resized_sim, metrics

########## Compyter vision Metrics ################
def calculate_fid(matrix1, matrix2):
    """
    Compute the Fréchet Inception Distance (FID) between two heatmaps.

    Parameters:
    matrix1 (numpy.ndarray): First matrix representing a 3C/Hi-C heatmap.
    matrix2 (numpy.ndarray): Second matrix representing a 3C/Hi-C heatmap.

    Returns:
    float: The FID score.
    """
    # Compute the mean and covariance of each matrix
    mu1, sigma1 = np.mean(matrix1, axis=0), np.cov(matrix1, rowvar=False)
    mu2, sigma2 = np.mean(matrix2, axis=0), np.cov(matrix2, rowvar=False)

    # Compute squared difference of means
    mean_diff = np.sum((mu1 - mu2) ** 2)

    # Compute sqrt of product of covariance matrices
    cov_sqrt = sqrtm(sigma1 @ sigma2)
    
    # Handle numerical issues
    if np.iscomplexobj(cov_sqrt):
        cov_sqrt = cov_sqrt.real

    # Compute FID score
    fid_score = mean_diff + np.trace(sigma1 + sigma2 - 2 * cov_sqrt)
    
    return fid_score

def compute_ssim(matrix1, matrix2):
    """
    Compute the Structural Similarity Index (SSIM) between two heatmaps.

    Parameters:
    matrix1 (numpy.ndarray): First heatmap.
    matrix2 (numpy.ndarray): Second heatmap.

    Returns:
    float: SSIM score (-1 to 1, higher is better).
    """
    return ssim(matrix1, matrix2, data_range=matrix2.max() - matrix2.min())

def compute_gmsd(matrix1, matrix2):
    """
    Compute the Gradient Magnitude Similarity Deviation (GMSD) between two heatmaps.

    Parameters:
    matrix1 (numpy.ndarray): First heatmap.
    matrix2 (numpy.ndarray): Second heatmap.

    Returns:
    float: GMSD score (lower means more similar).
    """
    # Compute gradients
    gx1 = cv2.Sobel(matrix1, cv2.CV_64F, 1, 0, ksize=3)
    gy1 = cv2.Sobel(matrix1, cv2.CV_64F, 0, 1, ksize=3)
    gx2 = cv2.Sobel(matrix2, cv2.CV_64F, 1, 0, ksize=3)
    gy2 = cv2.Sobel(matrix2, cv2.CV_64F, 0, 1, ksize=3)

    # Compute gradient magnitudes
    g1 = np.sqrt(gx1**2 + gy1**2)
    g2 = np.sqrt(gx2**2 + gy2**2)

    # Compute GMS
    C = 0.0026  # Small stability constant
    gms_map = (2 * g1 * g2 + C) / (g1**2 + g2**2 + C)
    
    # Compute standard deviation
    gmsd_value = np.std(gms_map)

    return gmsd_value

def compute_phase_correlation(matrix1, matrix2):
    """
    Compute the phase correlation similarity between two heatmaps.

    Parameters:
    matrix1 (numpy.ndarray): First heatmap.
    matrix2 (numpy.ndarray): Second heatmap.

    Returns:
    float: Peak phase correlation value (higher means more similar).
    """
    f1 = np.fft.fft2(matrix1)
    f2 = np.fft.fft2(matrix2)

    cross_power = (f1 * np.conj(f2)) / np.abs(f1 * np.conj(f2))  # Normalize
    correlation = np.fft.ifft2(cross_power).real

    return np.max(correlation)  # Peak similarity

########## Statistical Tests - Not the Best Methods #############

# 1. Mantel Test (Correlation Between Hi-C and Simulated Distance Matrices)
def mantel_test(distance_matrix, expected_dist_matrix, num_permutations=1000):
    """
    Compute the Mantel test between the simulated distance matrix and the expected one.

    Parameters:
    - distance_matrix (numpy.ndarray): Simulated 3D structure distance matrix.
    - expected_dist_matrix (numpy.ndarray): Expected distance matrix (derived from experimental Hi-C).
    - num_permutations (int): Number of permutations for significance testing.

    Returns:
    - observed_corr (float): Pearson correlation coefficient.
    - p_value (float): Significance of the correlation.
    """
    # Ensure diagonals are zero
    np.fill_diagonal(distance_matrix, 0)
    np.fill_diagonal(expected_dist_matrix, 0)

    # Make sure matrices are symmetric
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    expected_dist_matrix = (expected_dist_matrix + expected_dist_matrix.T) / 2

    # Convert to condensed form
    dist1 = ssd.squareform(distance_matrix)
    dist2 = ssd.squareform(expected_dist_matrix)

    # Compute observed Pearson correlation
    observed_corr = stats.pearsonr(dist1, dist2)[0]
    
    # Permutation test
    permuted_corrs = []
    for _ in range(num_permutations):
        permuted_indices = np.random.permutation(len(dist2))
        permuted_corr = stats.pearsonr(dist1, dist2[permuted_indices])[0]
        permuted_corrs.append(permuted_corr)

    # Compute p-value
    p_value = np.sum(np.abs(permuted_corrs) >= np.abs(observed_corr)) / num_permutations

    return observed_corr, p_value

# 2. Kullback-Leibler (KL) Divergence
def kl_divergence(distance_matrix, expected_dist_matrix):
    p_hist, _ = np.histogram(distance_matrix[np.triu_indices_from(distance_matrix, k=1)], bins=20, density=True)
    q_hist, _ = np.histogram(expected_dist_matrix[np.triu_indices_from(expected_dist_matrix, k=1)], bins=20, density=True)

    return stats.entropy(p_hist + 1e-10, q_hist + 1e-10)  # Avoid log(0)

# 3. Jensen-Shannon (JS) Divergence
def js_divergence(distance_matrix, expected_dist_matrix):
    p_hist, _ = np.histogram(distance_matrix[np.triu_indices_from(distance_matrix, k=1)], bins=20, density=True)
    q_hist, _ = np.histogram(expected_dist_matrix[np.triu_indices_from(expected_dist_matrix, k=1)], bins=20, density=True)

    return jensenshannon(p_hist, q_hist)

# 4. Earth Mover’s Distance (Wasserstein Distance)
def earth_movers_distance(distance_matrix, expected_dist_matrix):
    p_hist, _ = np.histogram(distance_matrix[np.triu_indices_from(distance_matrix, k=1)], bins=20, density=True)
    q_hist, _ = np.histogram(expected_dist_matrix[np.triu_indices_from(expected_dist_matrix, k=1)], bins=20, density=True)

    return wasserstein_distance(p_hist, q_hist)

# 5. Chi-Square Test
def chi_square_test(distance_matrix, expected_dist_matrix):
    p_hist, _ = np.histogram(distance_matrix[np.triu_indices_from(distance_matrix, k=1)], bins=20, density=True)
    q_hist, _ = np.histogram(expected_dist_matrix[np.triu_indices_from(expected_dist_matrix, k=1)], bins=20, density=True)

    # Normalize histograms to have the same sum
    p_hist = p_hist / np.sum(p_hist)
    q_hist = q_hist / np.sum(q_hist)

    return stats.chisquare(p_hist + 1e-10, q_hist + 1e-10)  # Avoid log(0)

# 6. Bhattacharyya Coefficient
def bhattacharyya_coefficient(distance_matrix, expected_dist_matrix):
    p_hist, _ = np.histogram(distance_matrix[np.triu_indices_from(distance_matrix, k=1)], bins=20, density=True)
    q_hist, _ = np.histogram(expected_dist_matrix[np.triu_indices_from(expected_dist_matrix, k=1)], bins=20, density=True)

    return np.sum(np.sqrt(p_hist * q_hist))

# 7. KS Test (Kolmogorov-Smirnov)
def ks_test_distance_distribution(distance_matrix, expected_dist_matrix):
    sim_distances = np.sort(distance_matrix[np.triu_indices_from(distance_matrix, k=1)])
    expected_distances = np.sort(expected_dist_matrix[np.triu_indices_from(expected_dist_matrix, k=1)])

    ks_stat, p_value = stats.ks_2samp(sim_distances, expected_distances)
    return ks_stat, p_value

def procrustes_test(matrix1, matrix2):
    """
    Perform Procrustes analysis to compare similarity between two matrices.

    Parameters:
    - matrix1, matrix2 (numpy.ndarray): Two matrices of the same shape.

    Returns:
    - float: Procrustes similarity score (lower means more similarity).
    """
    if matrix1.shape != matrix2.shape:
        raise ValueError("Matrices must have the same shape.")

    mtx1, mtx2, disparity = procrustes(matrix1, matrix2)
    return disparity

def rv_coefficient(matrix1, matrix2):
    """
    Compute the RV coefficient (a measure of similarity between two matrices).

    Parameters:
    - matrix1, matrix2 (numpy.ndarray): Two matrices of the same shape.

    Returns:
    - float: RV coefficient (higher means more similarity).
    """
    if matrix1.shape != matrix2.shape:
        raise ValueError("Matrices must have the same shape.")

    numerator = np.trace(matrix1.T @ matrix2 @ matrix2.T @ matrix1)
    denominator = np.sqrt(np.trace(matrix1.T @ matrix1 @ matrix1.T @ matrix1) *
                          np.trace(matrix2.T @ matrix2 @ matrix2.T @ matrix2))
    
    return numerator / denominator

def distance_correlation(X, Y):
    """
    Compute the distance correlation between two matrices.

    Parameters:
    - X, Y (numpy.ndarray): Two matrices.

    Returns:
    - float: Distance correlation (higher means more similarity).
    """
    def distance_matrix(A):
        return squareform(pdist(A, metric='euclidean'))

    def centering_matrix(D):
        n = D.shape[0]
        J = np.eye(n) - np.ones((n, n)) / n
        return J @ D @ J

    A = distance_matrix(X)
    B = distance_matrix(Y)

    A_centered = centering_matrix(A)
    B_centered = centering_matrix(B)

    dcov_XY = np.sqrt(np.sum(A_centered * B_centered) / (X.shape[0] ** 2))
    dcov_XX = np.sqrt(np.sum(A_centered * A_centered) / (X.shape[0] ** 2))
    dcov_YY = np.sqrt(np.sum(B_centered * B_centered) / (X.shape[0] ** 2))

    return dcov_XY / np.sqrt(dcov_XX * dcov_YY)

######## Stratified Distance Correlations #########

def compute_dscc(simulated_matrix, experimental_matrix, bin_size=5):
    """Compute distance-stratified correlation coefficients."""
    N = simulated_matrix.shape[0]
    max_distance = N // 2  # Avoid redundant comparisons
    genomic_distances = np.arange(1, max_distance, bin_size)
    
    pearson_correlations = []
    spearman_correlations = []

    for s in genomic_distances:
        indices = [(i, i + s) for i in range(N - s)]  # Select all pairs at distance s
        sim_distances = np.array([simulated_matrix[i, j] for i, j in indices])
        exp_distances = np.array([experimental_matrix[i, j] for i, j in indices])

        # Compute correlations if we have enough points
        if len(sim_distances) > 5:
            pearson_corr = stats.pearsonr(sim_distances, exp_distances)[0]
            spearman_corr = stats.spearmanr(sim_distances, exp_distances)[0]
        else:
            pearson_corr, spearman_corr = np.nan, np.nan  # Not enough data

        pearson_correlations.append(pearson_corr)
        spearman_correlations.append(spearman_corr)

    return genomic_distances, pearson_correlations, spearman_correlations

###### GenomeDISCO denoising technique ##########

def symmetric_normalization(A):
    """ Symmetric normalization: P = D^(-1/2) A D^(-1/2) """
    D = np.diag(1.0 / np.sqrt(A.sum(axis=1) + 1e-10))  # Avoid division by zero
    return D @ A @ D  # Equivalent to D^(-1/2) * A * D^(-1/2)

def column_stochastic_normalization(A):
    """ Column-stochastic normalization: P = A * D^(-1) where D is column sum matrix """
    col_sums = A.sum(axis=0, keepdims=True) + 1e-10  # Avoid division by zero
    return A / col_sums

def denoise_contact_map(A, t, method="symmetric"):
    """ Denoising using random walks with either symmetric or column normalization """
    if method == "symmetric":
        P = symmetric_normalization(A)
    elif method == "column":
        P = column_stochastic_normalization(A)
    else:
        raise ValueError("Unknown normalization method. Use 'symmetric' or 'column'.")

    P_t = np.linalg.matrix_power(P, t)  # Compute P^t
    return P_t

def compute_L1_difference(A1, A2):
    """ Compute the L1 distance between denoised contact maps A1 and A2. """
    # Count nonzero nodes in the original matrices
    N_nonzero = 0.5 * (np.count_nonzero(A1.sum(axis=1)) + np.count_nonzero(A2.sum(axis=1)))
    
    # Compute L1 distance
    L1_distance = np.sum(np.abs(A1 - A2)) / N_nonzero
    return L1_distance

def compute_concordance_score(A1, A2):
    """ Compute the concordance score R based on L1 difference. """
    L1_diff = compute_L1_difference(A1, A2)
    return 1 - L1_diff

####### HiC-Rep SCC metric ########

def mean_filter_smoothing(contact_map, h):
    """
    Applies 2D mean filter smoothing to a contact map.
    
    Args:
        contact_map (numpy array): The original contact map (n x n matrix).
        h (int): The smoothing window size parameter.
        
    Returns:
        numpy array: The smoothed contact map.
    """
    kernel_size = 1 + 2 * h  # Compute the smoothing window size
    smoothed_map = uniform_filter(contact_map, size=kernel_size, mode='nearest')
    return smoothed_map

def define_strata_from_matrices(matrix_shape, bin_size):
    """
    Generate strata (distance bins) based on matrix indices.
    
    Args:
        matrix_shape (tuple): Shape of the contact map matrix (n, n).
        bin_size (int): Bin size for stratification.

    Returns:
        numpy array: Stratum indices for each pair (i, j).
    """
    n = matrix_shape[0]
    strata = np.zeros((n, n), dtype=int)
    
    for i in range(n):
        for j in range(n):
            distance = abs(i - j)  # Compute genomic distance
            strata[i, j] = distance // bin_size  # Assign to bin
    
    return strata

def compute_scc(A1, A2, bin_size):
    """
    Compute the Stratum-Adjusted Correlation Coefficient (SCC) 
    using the CMH-based statistical approach.
    
    Args:
        A1 (numpy array): Contact map matrix for sample 1.
        A2 (numpy array): Contact map matrix for sample 2.
        bin_size (int): Bin size for stratification.
    
    Returns:
        float: SCC statistic (M^2), which should be in [-1,1].
    """
    assert A1.shape == A2.shape, "Matrices must have the same dimensions!"
    
    strata = define_strata_from_matrices(A1.shape, bin_size)
    unique_strata = np.unique(strata)
    
    N_k = []
    r_1k, r_2k = [], []
    rho_k = []
    var_Tk = []
    
    for k in unique_strata:
        mask = (strata == k)
        X_k = A1[mask]
        Y_k = A2[mask]
        N_k.append(len(X_k))

        if len(X_k) > 1:
            # Compute T_k
            T_k = np.sum(X_k * Y_k)
            
            # Compute E(T_k)
            E_Tk = (np.sum(X_k) * np.sum(Y_k)) / len(X_k)

            # Compute var(T_k)
            r1k = np.sum(X_k*Y_k)/len(X_k) - np.sum(X_k)*np.sum(Y_k)/len(X_k)**2
            r_1k.append(r1k)
            var_Xk = np.sum(X_k**2)/len(X_k)-(np.sum(X_k)/len(X_k))**2
            var_Yk = np.sum(Y_k**2)/len(Y_k)-(np.sum(Y_k)/len(Y_k))**2
            r2k = np.sqrt(var_Xk*var_Yk)
            r_2k.append(r2k)
            var_Tk_value = 1 / (len(X_k) - 1)*(np.sum(X_k**2)-np.sum(X_k)**2/len(X_k))*(np.sum(Y_k**2)-np.sum(Y_k)**2/len(Y_k)) if len(X_k) > 1 else 0
            var_Tk.append(var_Tk_value)
            rhok = r1k/r2k
        else:
            rhok = 0
            var_Tk.append(0)

        rho_k.append(rhok)

    N_k = np.array(N_k)
    r_1k, r_2k = np.array(r_1k), np.array(r_2k)
    rho_k = np.array(rho_k)
    var_Tk = np.array(var_Tk)

    # Compute rho_s as a weighted sum
    numerator = np.sum(N_k * r_2k * rho_k)
    denominator = np.sum(N_k * r_2k)
    rho_s = numerator / denominator if denominator != 0 else 0

    # Compute M^2 statistic
    M2 = (rho_s**2) * np.sum((N_k*r_2k)**2) / np.sum((N_k*r_2k)**2/(N_k-1)) if np.sum(N_k) != 0 else 0
    
    return M2, rho_s

def plot_heatmaps(matrix1, matrix2, matrix3, titles=('Data', 'Simulated', 'Random')):
    """
    Plot three heatmaps side-by-side.

    Parameters:
    matrix1, matrix2, matrix3 (numpy.ndarray): Heatmaps to plot.
    titles (tuple): Titles for each heatmap.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))  # Create 1 row, 3 columns
    
    heatmaps = [matrix1, matrix2, matrix3]
    
    for ax, heatmap, title in zip(axes, heatmaps, titles):
        im = ax.imshow(heatmap, cmap='Reds', aspect='auto',vmax=np.mean(heatmap)+np.std(heatmap))  # Heatmap with Viridis colormap
        ax.set_title(title)
        ax.axis('off')  # Remove axis labels
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)  # Add colorbar for each heatmap
    
    plt.tight_layout()
    plt.show()

def comparison_protocol(N_ensemble, ensemble_path, hic_path,label='multimm',diag_removal=False,N_diags=20):
    # Compute Simulated and Random heatmaps
    if label=='multimm':
        sim_heat = avg_heat_multimm(ensemble_path,N_ensemble)
    else:
        sim_heat = avg_heat_loopsage(ensemble_path,N_ensemble)
    sim_heat = (sim_heat-np.min(sim_heat))/(np.max(sim_heat)-np.min(sim_heat))
    _, rd_heat = generate_random_walks(N_ensemble,len(sim_heat))
    rd_heat = (rd_heat-np.min(rd_heat))/(np.max(rd_heat)-np.min(rd_heat))

    # Compute experimental heatmaps
    hic = hicstraw.HiCFile(hic_path)
    matrix_object = hic.getMatrixZoomData('1', '1', "observed", "NONE", "BP", 1000)
    exp_heat = matrix_object.getRecordsAsMatrix(178421513,179491193,178421513,179491193)
    exp_heat = (exp_heat-np.min(exp_heat))/(np.max(exp_heat)-np.min(exp_heat))

    # Compute power of power law
    p_s = average_signal_distance(exp_heat)
    s = np.arange(len(p_s))+1
    a, a_error = fit_power_law(s, p_s)
    print(f'Experimental power law fit: {a} =/- {a_error}')
    p_s = average_signal_distance(sim_heat)
    s = np.arange(len(p_s))+1
    a, a_error = fit_power_law(s, p_s)
    print(f'Simulated power law fit: {a} =/- {a_error}')
    p_s = average_signal_distance(rd_heat)
    s = np.arange(len(p_s))+1
    a, a_error = fit_power_law(s, p_s)
    print(f'Random power law fit: {a} =/- {a_error}')

    # Remove diagonals
    if diag_removal:
        rd_heat = remove_k_diagonals(rd_heat,N_diags)
        exp_heat = remove_k_diagonals(exp_heat,N_diags)
        sim_heat = remove_k_diagonals(sim_heat,N_diags)

    # Compute correlations
    print('\n---- Correlations -----')
    exp_new, rd_new, metrics = compare_heatmaps(exp_heat, rd_heat)
    print('Correlations of random model.',metrics)
    exp_new, sim_new, metrics = compare_heatmaps(exp_heat, sim_heat)
    print('Correlations of simulation model.',metrics)
    plot_heatmaps(exp_new, sim_new, rd_new)

    # Computer Vision metrics
    print('\n---- Computer Vision Metrics -----')
    fid_value = calculate_fid(exp_new, rd_new)
    print("FID Score (random):", fid_value)
    fid_value = calculate_fid(exp_new, sim_new)
    print("FID Score (simulation):", fid_value)

    ssim = compute_ssim(exp_new, rd_new)
    print("\nSSIM Score (random):", ssim)
    ssim = compute_ssim(exp_new, sim_new)
    print("SSIM Score (simulation):", ssim)

    gmsd = compute_gmsd(exp_new, rd_new)
    print("\nGMSD Score (random):", gmsd)
    gmsd = compute_gmsd(exp_new, sim_new)
    print("GMSD Score (simulation):", gmsd)

    phase_corr_score = compute_phase_correlation(exp_new, rd_new)
    print("\nPhase Correlation Score (random):", phase_corr_score)
    phase_corr_score = compute_phase_correlation(exp_new, sim_new)
    print("Phase Correlation Score (simulation):", phase_corr_score)

    # Other statistical tests
    ## Tests with simulation
    mantel_corr, mantel_p = mantel_test(sim_new, exp_new)
    kl_div = kl_divergence(sim_new, exp_new)
    js_div = js_divergence(sim_new, exp_new)
    emd = earth_movers_distance(sim_new, exp_new)
    chi2_stat, chi2_p = chi_square_test(sim_new, exp_new)
    b_coeff = bhattacharyya_coefficient(sim_new, exp_new)
    ks_stat, ks_p_value = ks_test_distance_distribution(sim_new, exp_new)
    proc_stat = procrustes_test(sim_new, exp_new)
    dist_corr = distance_correlation(sim_new, exp_new)
    rv_coff = rv_coefficient(sim_new, exp_new)

    print('\n---Tests with Simulated Model---')
    print(f"Mantel Test: Correlation={mantel_corr:.3f}, p-value={mantel_p:.3f}")
    print(f"KL Divergence: {kl_div:.3f}")
    print(f"JS Divergence: {js_div:.3f}")
    print(f"Earth Mover’s Distance: {emd:.3f}")
    print(f"Chi-Square Test: Statistic={chi2_stat:.3f}, p-value={chi2_p:.3f}")
    print(f"Bhattacharyya Coefficient: {b_coeff:.3f}")
    print(f"KS Test Statistic: {ks_stat:.3f}, p-value={ks_p_value:.3f}")
    print(f"Procrustes stat: {proc_stat}")
    print(f"Distance Correlation: {dist_corr}")
    print(f"RV coefficient: {rv_coff}")

    ## Tests with random models
    mantel_corr, mantel_p = mantel_test(rd_new, exp_new)
    kl_div = kl_divergence(rd_new, exp_new)
    js_div = js_divergence(rd_new, exp_new)
    emd = earth_movers_distance(rd_new, exp_new)
    chi2_stat, chi2_p = chi_square_test(rd_new, exp_new)
    b_coeff = bhattacharyya_coefficient(rd_new, exp_new)
    ks_stat, ks_p_value = ks_test_distance_distribution(rd_new, exp_new)
    proc_stat = procrustes_test(rd_new, exp_new)
    dist_corr = distance_correlation(rd_new, exp_new)
    rv_coff = rv_coefficient(rd_new, exp_new)

    print('\n----Tests with Random Model----')
    print(f"Mantel Test: Correlation={mantel_corr:.3f}, p-value={mantel_p:.3f}")
    print(f"KL Divergence: {kl_div:.3f}")
    print(f"JS Divergence: {js_div:.3f}")
    print(f"Earth Mover’s Distance: {emd:.3f}")
    print(f"Chi-Square Test: Statistic={chi2_stat:.3f}, p-value={chi2_p:.3f}")
    print(f"Bhattacharyya Coefficient: {b_coeff:.3f}")
    print(f"KS Test Statistic: {ks_stat:.3f}, p-value={ks_p_value:.3f}")
    print(f"Procrustes stat: {proc_stat}")
    print(f"Distance Correlation: {dist_corr}")
    print(f"RV coefficient: {rv_coff}")

    # Stratified Correlation
    genomic_distances, pearson_corrs, spearman_corrs = compute_dscc(sim_new, exp_new)
    genomic_distances, pearson_corrs_rd, spearman_corrs_rd = compute_dscc(rd_new, exp_new)
    plt.figure(figsize=(15, 7))
    plt.plot(genomic_distances, pearson_corrs, label="Pearson (simulation)", marker="o" , color='orange')
    plt.plot(genomic_distances, spearman_corrs, label="Spearman (simulation)", marker="s", color='red')
    plt.plot(genomic_distances, pearson_corrs_rd, label="Pearson (random walk)", marker="o" , color='cyan')
    plt.plot(genomic_distances, spearman_corrs_rd, label="Spearman (ramdom walk)", marker="s", color='blue')
    plt.xlabel("Genomic Distance (bins)",fontsize=18)
    plt.ylabel("Correlation Coefficient",fontsize=18)
    plt.legend()
    plt.title("Distance-Stratified Correlation Coefficients",fontsize=23)
    plt.grid()
    plt.show()
    
    # GenomeDISCO approach
    t=3
    denoised_sim = denoise_contact_map(sim_new, t)
    denoised_exp = denoise_contact_map(exp_new, t)
    denoised_rd = denoise_contact_map(rd_new, t)
    concordance_score = compute_concordance_score(denoised_sim, denoised_exp)
    print(f"\nConcordance Score (simulation versus experiment): {concordance_score}")
    concordance_score = compute_concordance_score(denoised_rd, denoised_exp)
    print(f"Concordance Score (random versus experiment): {concordance_score}")

    # HiC-Rep approach
    _, scc_sim = compute_scc(sim_new, exp_new, 200)
    _, scc_rd = compute_scc(rd_new, exp_new, 200)
    print('\nSCC with simulated data:',scc_sim)
    print('SCC with random data:',scc_rd)

def compare_heats(matrix1,matrix2,diag_removal=False,N_diags=40):
    # Compute power of power law
    p_s = average_signal_distance(matrix1)
    s = np.arange(len(p_s))+1
    a, a_error = fit_power_law(s, p_s)
    print(f'Experimental power law fit: {a} =/- {a_error}')
    p_s = average_signal_distance(matrix2)
    s = np.arange(len(p_s))+1
    a, a_error = fit_power_law(s, p_s)
    print(f'Simulated power law fit: {a} =/- {a_error}')

    # Remove diagonals
    if diag_removal:
        matrix1 = remove_k_diagonals(matrix1,N_diags)
        matrix2 = remove_k_diagonals(matrix2,N_diags)

    # Compute correlations
    print('\n---- Correlations -----')
    mat1_new, mat2_new, metrics = compare_heatmaps(matrix1,matrix2)
    print('Correlations of random model.',metrics)

    # Computer Vision metrics
    print('\n---- Computer Vision Metrics -----')
    fid_value = calculate_fid(mat1_new, mat2_new)
    print("FID Score (random):", fid_value)

    ssim = compute_ssim(mat1_new, mat2_new)
    print("\nSSIM Score (random):", ssim)

    gmsd = compute_gmsd(mat1_new, mat2_new)
    print("\nGMSD Score (random):", gmsd)

    phase_corr_score = compute_phase_correlation(mat1_new, mat2_new)
    print("\nPhase Correlation Score (random):", phase_corr_score)

    # Other statistical tests
    ## Tests with simulation
    mantel_corr, mantel_p = mantel_test(mat1_new, mat2_new)
    kl_div = kl_divergence(mat1_new, mat2_new)
    js_div = js_divergence(mat1_new, mat2_new)
    emd = earth_movers_distance(mat1_new, mat2_new)
    chi2_stat, chi2_p = chi_square_test(mat1_new, mat2_new)
    b_coeff = bhattacharyya_coefficient(mat1_new, mat2_new)
    ks_stat, ks_p_value = ks_test_distance_distribution(mat1_new, mat2_new)
    proc_stat = procrustes_test(mat1_new, mat2_new)
    dist_corr = distance_correlation(mat1_new, mat2_new)
    rv_coff = rv_coefficient(mat1_new, mat2_new)

    print('\n---Tests with Simulated Model---')
    print(f"Mantel Test: Correlation={mantel_corr:.3f}, p-value={mantel_p:.3f}")
    print(f"KL Divergence: {kl_div:.3f}")
    print(f"JS Divergence: {js_div:.3f}")
    print(f"Earth Mover’s Distance: {emd:.3f}")
    print(f"Chi-Square Test: Statistic={chi2_stat:.3f}, p-value={chi2_p:.3f}")
    print(f"Bhattacharyya Coefficient: {b_coeff:.3f}")
    print(f"KS Test Statistic: {ks_stat:.3f}, p-value={ks_p_value:.3f}")
    print(f"Procrustes stat: {proc_stat}")
    print(f"Distance Correlation: {dist_corr}")
    print(f"RV coefficient: {rv_coff}")

    # Stratified Correlation
    genomic_distances, pearson_corrs, spearman_corrs = compute_dscc(mat1_new, mat2_new)
    plt.figure(figsize=(15, 7))
    plt.plot(genomic_distances, pearson_corrs, label="Pearson", marker="o" , color='orange')
    plt.plot(genomic_distances, spearman_corrs, label="Spearman", marker="s", color='red')
    plt.xlabel("Genomic Distance (bins)",fontsize=18)
    plt.ylabel("Correlation Coefficient",fontsize=18)
    plt.legend()
    plt.title("Distance-Stratified Correlation Coefficients",fontsize=23)
    plt.grid()
    plt.show()
    
    # GenomeDISCO approach
    t=3
    denoised1 = denoise_contact_map(mat1_new, t)
    denoised2 = denoise_contact_map(mat2_new, t)
    concordance_score = compute_concordance_score(denoised1, denoised2)
    print(f"\nConcordance Score (simulation versus experiment): {concordance_score}")

    # HiC-Rep approach
    _, scc_sim = compute_scc(mat1_new, mat2_new, 200)
    print('\nSCC with simulated data:',scc_sim)

def main():
    N_ensemble = 950
    label = 'loopsage'
    ensemble_path = '/home/skorsak/Data/hackathon/LoopSage_ensemble_Annealing_two_fams'
    hic_path = '/home/skorsak/Data/4DNucleome/4DNFIC3JD6O2/4DNFIC3JD6O2.hic'
    comparison_protocol(N_ensemble, ensemble_path, hic_path, label,diag_removal=False,N_diags=50)

if __name__=='__main__':
    main()