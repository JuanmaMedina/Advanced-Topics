import numpy as np
from scipy.stats import chisquare

# Load bases and their QSs
bases = open("/home/juanma/TOPICS/NGS_EM/Exercise1/bases.txt", 'r').readlines()
qs = open("/home/juanma/TOPICS/NGS_EM/Exercise1/qualities.txt", 'r').readlines()

# Iterate over bases and qs, parse them and introduce the qs of each read into a list
for i in np.argsort(bases):
    bases[i] = bases[i].replace(" ", "").rstrip()
    qs[i] = qs[i].rstrip().split(" ")

# Iterate over chromosomic positions
for i in np.argsort(bases):

    # Store base frequencies in empty array: A, C, G, T
    base_freq = np.zeros(shape = (1, 4))

    # Initial base-frequencies random guess (A-T, C-G coupled)
    initial_freq = np.array([0.3, 0.2, 0.2, 0.3])

    # Sanity check
    if sum(initial_freq) != 1:
        raise ValueError('Base-frequencies must sum up to one!')

    # Store genotype LHs in empty lists
    gen_A, gen_C, gen_G, gen_T = [], [], [], []

    # Iterate over bases
    for j in range(0, len(bases[i])):

        # LH estimation using the GATK model:

        # If the observed base (X) does not match the real base (Z), add the error probability (Z = E / 3) as a float
        # to the corresponding genotype LHs lists
        # If X matches Z, add the probability (Z = 1 - E) to the genotype LH list of that base

        if bases[i][j] == "A":
            gen_C.append(float(qs[i][j]) / 3), gen_G.append(float(qs[i][j]) / 3), gen_T.append(float(qs[i][j]) / 3)
            gen_A.append(1 - float(qs[i][j]))

        elif bases[i][j] == "C":
            gen_A.append(float(qs[i][j]) / 3), gen_G.append(float(qs[i][j]) / 3), gen_T.append(float(qs[i][j]) / 3)
            gen_C.append(1 - float(qs[i][j]))

        elif bases[i][j] == "G":
            gen_A.append(float(qs[i][j]) / 3), gen_C.append(float(qs[i][j]) / 3), gen_T.append(float(qs[i][j]) / 3)
            gen_G.append(1 - float(qs[i][j]))

        else:
            gen_A.append(float(qs[i][j]) / 3), gen_C.append(float(qs[i][j]) / 3), gen_G.append(float(qs[i][j]) / 3)
            gen_T.append(1 - float(qs[i][j]))


    # EM-convergence condition
    convergence = abs(base_freq - initial_freq)

    # EM iterates until the inter-step distance is less than 0.001
    while np.min(convergence) > 0.001:

        # Base-frequencies update
        base_freq = initial_freq

        # Store marginal and overall expectations in empty lists
        A_exp, C_exp, G_exp, T_exp, total_exp = [], [], [], [], []

        # E-step: iterate over bases and multiply by their genotype LHs
        for z in range(0, len(bases[i])):

            # Marginal expectations (numerator: prior x LH)
            A_exp.append(base_freq[0] * gen_A[z])
            C_exp.append(base_freq[1] * gen_C[z])
            G_exp.append(base_freq[2] * gen_G[z])
            T_exp.append(base_freq[3] * gen_T[z])

            # Overall expectation (denominator)
            total_exp.append(base_freq[0] * gen_A[z] + base_freq[1] * gen_C[z] +
                             base_freq[2] * gen_G[z] + base_freq[3] * gen_T[z])

        # M-step: mean of the marginal expectations divided by the total expectation
        initial_freq = np.array([np.mean(np.array(A_exp)/np.array(total_exp)),
                                 np.mean(np.array(C_exp)/np.array(total_exp)),
                                 np.mean(np.array(G_exp)/np.array(total_exp)),
                                 np.mean(np.array(T_exp)/np.array(total_exp))])

        # Convergence update
        convergence = abs(base_freq - initial_freq)


# Allele frequencies
print("The allele frequencies are:")
print("A: %s" % base_freq[0])
print("C: %s" % base_freq[1])
print("G: %s" % base_freq[2])
print("T: %s" % base_freq[3])


