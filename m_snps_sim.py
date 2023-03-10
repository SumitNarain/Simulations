from random import random

import matplotlib.pyplot as plt
import numpy as np
import scipy as sc
import math



def m_snps_sim(rs_ids, sim_correlation):

    #Obtains the first row of the sim correlation and sets sample size (n)
    M = sim_correlation.shape[0]
    n = 100000

    maf_v = np.random.uniform(0, 1, M)
    genotype = np.random.binomial(2, maf_v, (n, M)).T
    # samples = list()
    # for j in range(M):
    #     aux = np.random.binomial(2, maf_v[j], n)
    #     mean = np.mean(aux)
    #     std = np.std(aux)
    #     # samples.append(aux)
    #     samples.append((aux - mean) / std)
    # genotype = np.asarray(samples)

    L = np.linalg.cholesky(sim_correlation)
    corr_x = np.matmul(L, genotype)
    for j in range(M):
        corr_x[j] = (corr_x[j] - corr_x[j].mean()) / corr_x[j].std()

    # assigment of heritability and sigma_e assuming var[y] = 1
    h_2 = 0.5
    p = 0.1
    M_c = int(math.ceil(M * p))
    sigma_e = 1 - h_2
    true_beta = np.concatenate((np.random.normal(0, np.sqrt(h_2/M_c), M_c), np.random.normal(0, 0.0000000001, M-M_c)))
    Xb = np.matmul(corr_x.T, true_beta)
    noise = np.random.normal(0, np.sqrt(sigma_e), n)
    y = Xb + noise
    y = (y - y.mean())

    h_2_calc = Xb.var() / (Xb.var() + noise.var())
    heritability = sum(true_beta ** 2)

    noise_var_calc = noise.var() / (Xb.var() + noise.var())

    marginal_beta = np.zeros(len(maf_v))
    marginal_var = np.zeros(len(maf_v))
    for j in range(M):
        marginal_beta[j] = np.matmul(corr_x[j].reshape((1,n)), y.reshape((n,1))) / n
        has_to_be_n = np.matmul(corr_x[j].T, corr_x[j])
        marginal_var[j] = sum((y - (corr_x[j]*marginal_beta[j]))**2) / ((n - 1) * has_to_be_n)
    marginal_se = np.sqrt(marginal_var)

    marginal_z = marginal_beta / marginal_se

    plt.figure(1)
    plt.scatter(true_beta, marginal_beta)
    plt.title("true betas vs marginal betas")
    plt.figure(2)
    plt.hist(true_beta)
    plt.title("true betas histogram")
    plt.figure(3)
    plt.hist(marginal_beta)
    plt.title("marginal betas histogram")
    plt.show()

    marginal_p = np.zeros(M)
    for j in range(len(maf_v)):
        marginal_p[j] = sc.stats.norm.sf(abs(marginal_z[j])) * 2

    # rs_ids = ['rs{}'.format(val) for val in range(M)]
    file_name = 'sim_M_{}'.format(M)
    file = open(file_name, 'w')
    file.write("MarkerName\tEffect\tStdErr\tGC.Zscore\tGC.Pvalue\tWeight\tAllele1\tAllele2\n")
    for j in range(len(maf_v)):
        file.write(
            "{}\t{}\t{}\t{}\t{}\t{}\tC\tG\n".format(rs_ids[j], marginal_beta[j], marginal_se[j], marginal_z[j], marginal_p[j],
                                              n))

    print("Heritability: {}".format(heritability))
    file.close()
    return rs_ids, marginal_beta, marginal_se, marginal_z, marginal_p, n, M, true_beta, h_2_calc, sim_correlation


if __name__ == "__main__":
    #generating random variant names and creating a region correlation matrix of [100,100]
    snp_rsid = ["rs" + str(x) for x in range(10000)]
    region_correlation = np.diag(np.ones(10000))
    rs_ids, marginal_beta, marginal_se, marginal_z, marginal_p, n, M, true_beta, heritability, correlation = m_snps_sim(
        snp_rsid, region_correlation)




