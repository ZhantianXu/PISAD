#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import scipy
from scipy.stats import nbinom
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import math
from pathlib import Path
import os

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="GenomeScope-like model fitting for k-mer frequency histograms.")
    parser.add_argument("-i", "--input", required=True, help="Path to the k-mer histogram file.")
    parser.add_argument("-m", "--model", type=float, required=True, help="Heterozygosity ratio (0 or 1).")

    # Optional parameters
    parser.add_argument("-o", "--output", help="Path to save the output plot.")
    parser.add_argument("-est", "--est_kmercov", type=int, help="Estimated k-mer coverage (optional).")
    parser.add_argument("-cutoff", "--cutoff_threshold", type=float, help="Coverage threshold (0.95)(optional).")
    parser.add_argument("-het", "--het_ratio", type=float, help="Initial heterozygosity(0/0.12) (optional).")
    parser.add_argument("-rho", "--rho", type=float, help="Initial rho value (0.2)(optional).")
    parser.add_argument("-k", "--kmer", type=int, help="kmer(21)(optional)")

    return parser.parse_args()


def load_kmer_histogram_all(hist_file):
    """Load k-mer frequency histogram file."""
    hist_df = pd.read_csv(hist_file, sep='\t', header=None, names=['coverage', 'frequency'], skiprows=0, nrows=1000)
    coverage = hist_df['coverage'].values
    frequency = hist_df['frequency'].values
    return coverage, frequency


def negative_binomial_pmf(x, mu, size):
    """Negative binomial probability mass function (PMF)."""
    return nbinom.pmf(x, n=size, p=size / (size + mu))


def genomescope_diploid_model_nb(x, r, d, lambda_val, rho, G, kmer_length=21):
    """GenomeScope diploid mixture model function."""
    alpha_val = 2 * (1 - d) * (1 - (1 - r) ** kmer_length) + 2 * d * (1 - (1 - r) ** kmer_length) ** 2 + 2 * d * (
                (1 - r) ** kmer_length) * (1 - (1 - r) ** kmer_length)
    beta_val = (1 - d) * ((1 - r) ** kmer_length) + d * (1 - (1 - r) ** kmer_length) ** 2
    term1 = alpha_val * negative_binomial_pmf(x, lambda_val, lambda_val / rho)
    term2 = beta_val * negative_binomial_pmf(x, 2 * lambda_val, 2 * lambda_val / rho)
    return G * (term1 + term2)

def genomescope_diploid_model_nb1(x, r, d, lambda_val, rho, G, kmer_length=21):
    """GenomeScope diploid mixture model function (revised)."""
    alpha_val = 2 * (1 - d) * (1 - (1 - r) ** kmer_length) + 2 * d * (1 - (1 - r) ** kmer_length) ** 2 + 2 * d * (
                (1 - r) ** kmer_length) * (1 - (1 - r) ** kmer_length)
    term1 = alpha_val * negative_binomial_pmf(x, lambda_val, lambda_val / rho)
    return G * term1


def genomescope_diploid_model_nb2(x, r, d, lambda_val, rho, G, kmer_length=21):
    """GenomeScope diploid mixture model function (revised)."""
    beta_val = (1 - d) * ((1 - r) ** kmer_length) + d * (1 - (1 - r) ** kmer_length) ** 2
    term2 = beta_val * negative_binomial_pmf(x, 2 * lambda_val, 2 * lambda_val / rho)
    return G * term2


def estimate_initial_params(coverage_data, frequency_data, model,threshold):
    """Estimate estKmerCov and estGenomeSize."""
    cumulative_sum = np.cumsum(frequency_data)
    total_sum = cumulative_sum[-1]

    threshold_index = np.searchsorted(cumulative_sum, threshold * total_sum)
    max_coverage_95 = coverage_data[min(threshold_index, len(frequency_data) - 1)]

    peak_indices = []
    for i in range(1, len(frequency_data) - 1):
        if frequency_data[i - 1] < frequency_data[i] and frequency_data[i] > frequency_data[i + 1]:
            if cumulative_sum[i] / total_sum <= threshold:
                peak_indices.append(i)

    if len(peak_indices) == 1:
        est_kmercov = int(coverage_data[peak_indices[0]] / 2) if model == 0 else int(coverage_data[peak_indices[0]])
    elif len(peak_indices) == 2:
        peak1, peak2 = coverage_data[peak_indices[0]], coverage_data[peak_indices[1]]
        ratio = peak2 / peak1
        est_kmercov = peak1 if 1.8 <= ratio <= 2.2 else (int(peak1 / 2) if model == 0 else int(peak1))
    else:
        est_kmercov = 1 if model == 0 else 2
        print(f"The coverage is too low, the selection of hybrid peak may be wrong !!")

    total_kmers = np.sum(frequency_data)
    est_genomesize = total_kmers / est_kmercov

    return max_coverage_95, est_kmercov, est_genomesize


if __name__ == "__main__":
    # Parse arguments
    args = parse_arguments()
    input_file = args.input
    file_name = Path(input_file).stem
    model = args.model

    # Heterozygosity bounds
    min_het, max_het = (0, 0.012) if model == 0 else (0.012, 1)

    # Load k-mer data
    coverage_all, frequency_all = load_kmer_histogram_all(input_file)

    threshold = args.cutoff_threshold if args.cutoff_threshold is not None else 0.95

    # Estimate initial parameters
    max_coverage_95, default_est_kmercov, est_genomesize = estimate_initial_params(coverage_all, frequency_all, model,threshold)

    # Use provided parameters or defaults
    est_kmercov = args.est_kmercov if args.est_kmercov is not None else default_est_kmercov
    initial_het = args.het_ratio if args.het_ratio is not None else (min_het + 1e-3)
    initial_rho = args.rho if args.rho is not None else 0.2
    k = args.kmer if args.kmer is not None else 21

    print(f"Using parameters: threshold={threshold},cut_off={max_coverage_95},est_kmercov={est_kmercov}, initial_het={initial_het}, initial_rho={initial_rho}")

    # Initial guess for curve fitting
    initial_params_guess = [initial_het, 0, est_kmercov, initial_rho, est_genomesize]

    # Fit the model
    optimal_params, _ = curve_fit(
        genomescope_diploid_model_nb,
        coverage_all[math.floor(est_kmercov / 2):],
        frequency_all[math.floor(est_kmercov / 2):],
        p0=initial_params_guess,
        bounds=([min_het, 0, 0, 0.15, 0], [max_het, 1, np.inf, 0.65, np.inf])
    )
    r_opt, d_opt, lambda_opt, rho_opt, G_opt = optimal_params

    print("\nFitting successful!")
    print(f"r (heterozygosity): {r_opt:.4f}")
    print(f"d (repeat fraction): {d_opt:.4f}")
    print(f"lambda (haploid coverage): {lambda_opt:.2f}")
    print(f"G (genome scaling factor): {G_opt:.2e}")
    print(f"rho: {rho_opt:.2e}")

    # Plot results
    x_plot = np.linspace(1, int(lambda_opt * 5), int(lambda_opt * 5))
    y_predicted1 = genomescope_diploid_model_nb1(x_plot, *optimal_params)
    y_predicted2 = genomescope_diploid_model_nb2(x_plot, *optimal_params)
    y_predicted = genomescope_diploid_model_nb(x_plot, *optimal_params)

    # Compute residuals in the range (0, est_kmercov)
    mask_diff = (coverage_all >= 0) & (coverage_all <= est_kmercov)
    x_diff = coverage_all[mask_diff]
    observed_diff = frequency_all[mask_diff]
    predicted_diff = genomescope_diploid_model_nb(x_diff, *optimal_params)
    residual = observed_diff - predicted_diff
    residual = np.clip(residual, 0, None)  # Replace all negative values with 0

    total_error_kmers = np.sum(residual * x_diff)
    total_kmers = np.sum(coverage_all * frequency_all)
    error_rate = 1 - (1 - (total_error_kmers / total_kmers)) ** (1 / k)
    genome_length = (total_kmers - total_error_kmers) / (lambda_opt * 2)

    length_diff = len(y_predicted) - len(residual)
    residual_padded = np.pad(residual, (0, length_diff), mode='constant', constant_values=0)
    r2 = r2_score(frequency_all[int(lambda_opt):int(lambda_opt * 5)], y_predicted[int(lambda_opt):])

    left = int(np.floor(-est_kmercov/2 + est_kmercov))
    if(left<2):
        left=2
    right = int(np.ceil(est_kmercov/2 + est_kmercov))

    print(f"\nleft,right: {left} - {right}")
    print(f"error_rate: {error_rate:.4f}")
    print(f"genome_length: {genome_length:.2e}")

    plt.figure(figsize=(10, 6))
    plt.bar(coverage_all[:int(lambda_opt * 5)], frequency_all[:int(lambda_opt * 5)], label='observed',
            width=0.9, alpha=0.7)
    plt.axvline(x=max_coverage_95, color='r', linestyle='--', linewidth=2, label='cut-off threshold')
    plt.axvline(x=lambda_opt, color='black', linestyle='--', linewidth=1, label='peaks')
    plt.axvline(x=lambda_opt*2, color='black', linestyle='--', linewidth=1)
    plt.plot(x_plot, y_predicted1, '-', label='heterozygous', linewidth=2, color='orange')
    plt.plot(x_plot, y_predicted2, '-', label='homozygous', linewidth=2, color='yellow')
    if len(x_diff) > 1:
        plt.plot(x_diff, residual, label='erroneous', linewidth=2, color='green')
    else:
        plt.vlines(x=x_diff[0], ymin=0, ymax=residual[0], colors='green', linestyles='--', linewidth=2,
                   label='Residual')
    plt.plot(x_plot, y_predicted+residual_padded, '-', label='full model', linewidth=3, color='black')

    title = (f'Diploid Model Fitting\n'
             '\n'
             f'len: {genome_length:.2e}    '
             f'r: {r_opt*100:.3f}%    '
             f'err: {error_rate*100:.3f}%    '
             f'kcov: {lambda_opt:.1f}    '
             f'R^2: {r2:.2f}    '
             f'range: {left}-{right}'
             )

    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend()
    plt.ylim(0, max(frequency_all[est_kmercov*2]*2,frequency_all[est_kmercov*1]*2))
    plt.gca().set_facecolor('lightgray')
    plt.axvspan(left, right, color='gray', alpha=0.3, label='Marked Area')
    output_file = f"{file_name}_{left}_{right}.png"
    if args.output:
        output_file = os.path.join(args.output, output_file)
    plt.savefig(output_file)
    print(f"Graph saved to {output_file}")

