from scipy.linalg import expm, logm


# Loading pedigree data
def load_pedigree(path):
    return np.loadtxt(path, delimiter='\t', skiprows=1)


# Calculating initial genotype probabilities for homozygous and heterozygous
def initial_probabilities(base_probabilities, prop_het):
    prob_aa, prob_cc, prob_tt, prob_gg = base_probabilities
    p0aa = prob_aa * (1 - prop_het)
    p0cc = prob_cc * (1 - prop_het)
    p0tt = prob_tt * (1 - prop_het)
    p0gg = prob_gg * (1 - prop_het)

    heterozygous_probs = prop_het / 12
    p0ac = heterozygous_probs
    p0at = heterozygous_probs
    p0ag = heterozygous_probs
    p0ca = heterozygous_probs
    p0ct = heterozygous_probs
    p0cg = heterozygous_probs
    p0ta = heterozygous_probs
    p0tc = heterozygous_probs
    p0tg = heterozygous_probs
    p0ga = heterozygous_probs
    p0gc = heterozygous_probs
    p0gt = heterozygous_probs

    initial_state_probs = np.array(
        [p0aa, p0cc, p0tt, p0gg, p0ac, p0at, p0ag, p0ca, p0ct, p0cg, p0ta, p0tc, p0tg, p0ga, p0gc, p0gt])

    initial_state_probs /= initial_state_probs.sum()

    return initial_state_probs


# Constructing the 16x16 transition matrix
def transition_matrix_construction(g):
    # homozygous → homozygous (4x4)
    mat_T1 = np.array([
        [(1 - g) ** 2, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9],
        [g ** 2 / 9, (1 - g) ** 2, g ** 2 / 9, g ** 2 / 9],
        [g ** 2 / 9, g ** 2 / 9, (1 - g) ** 2, g ** 2 / 9],
        [g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, (1 - g) ** 2]
    ])

    # homozygous → heterozygous (4x12)
    mat_T2 = np.array([
        [1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g,
         g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9,
         g ** 2 / 9],
        [1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g,
         1 / 3 * (1 - g) * g, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9],
        [g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9,
         1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g],
        [g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g,
         g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g]
    ])

    # heterozygous → homozygous (12x4)
    mat_T3 = np.array([
        [1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9],
        [1 / 3 * (1 - g) * g, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9],
        [1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g],
        [1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9],
        [g ** 2 / 9, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, g ** 2 / 9],
        [g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, 1 / 3 * (1 - g) * g],
        [1 / 3 * (1 - g) * g, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9],
        [g ** 2 / 9, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, g ** 2 / 9],
        [g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g],
        [1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g],
        [g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, 1 / 3 * (1 - g) * g],
        [g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g]
    ])

    # heterozygous → heterozygous (12x12)
    mat_T4 = np.array([
        [(1 - g) ** 2, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9,
         g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9],
        [1 / 3 * (1 - g) * g, (1 - g) ** 2, 1 / 3 * (1 - g) * g, g ** 2 / 9,
         1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9,
         1 / 3 * (1 - g) * g],
        [1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, (1 - g) ** 2, g ** 2 / 9, g ** 2 / 9,
         1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9],
        [g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, (1 - g) ** 2, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g,
         1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9],
        [g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, 1 / 3 * (1 - g) * g, (1 - g) ** 2,
         1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g],
        [g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g,
         (1 - g) ** 2, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9],
        [g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9,
         (1 - g) ** 2, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9],
        [1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9,
         1 / 3 * (1 - g) * g, (1 - g) ** 2, 1 / 3 * (1 - g) * g, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9],
        [g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g,
         1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, (1 - g) ** 2, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9],
        [g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9,
         1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, (1 - g) ** 2, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g],
        [1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, g ** 2 / 9,
         1 / 3 * (1 - g) * g, g ** 2 / 9, 1 / 3 * (1 - g) * g, (1 - g) ** 2, 1 / 3 * (1 - g) * g],
        [g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, g ** 2 / 9,
         g ** 2 / 9, g ** 2 / 9, g ** 2 / 9, 1 / 3 * (1 - g) * g, 1 / 3 * (1 - g) * g, (1 - g) ** 2]
    ])

    # Constructing the whole matrix G (16x16)
    mat_G = np.zeros((16, 16))
    mat_G[0:4, 0:4] = mat_T1
    mat_G[0:4, 4:] = mat_T2
    mat_G[4:, 0:4] = mat_T3
    mat_G[4:, 4:] = mat_T4

    return mat_G


# Divergence effects matrix: each entry tells how different genotypes are from each other
def get_divergence_effects():
    return np.array([
        [0, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 1, 1, 0.5, 1, 1, 0.5, 1, 1],
        [1, 0, 1, 1, 0.5, 1, 1, 0.5, 1, 1, 1, 0.5, 1, 1, 0.5, 1],
        [1, 1, 0, 1, 1, 0.5, 1, 1, 0.5, 1, 1, 1, 0.5, 1, 1, 1],
        [1, 1, 1, 0, 1, 1, 0.5, 1, 1, 0.5, 1, 1, 1, 0.5, 1, 1],
        [0.5, 0.5, 1, 1, 0, 1, 1, 1, 0.5, 1, 1, 1, 0.5, 0.5, 1, 0.5],
        [0.5, 1, 0.5, 1, 1, 0, 1, 1, 1, 1, 0.5, 1, 1, 1, 0.5, 0.5],
        [0.5, 1, 1, 0.5, 1, 1, 0, 1, 1, 1, 1, 0.5, 1, 1, 0.5, 1],
        [0.5, 0.5, 1, 1, 1, 1, 1, 0, 1, 0.5, 1, 1, 0.5, 1, 1, 0.5],
        [1, 1, 0.5, 1, 0.5, 1, 1, 1, 0, 1, 1, 1, 0.5, 0.5, 0.5, 1],
        [1, 1, 1, 0.5, 1, 1, 1, 0.5, 1, 0, 0.5, 1, 1, 1, 1, 1],
        [0.5, 1, 1, 1, 1, 0.5, 1, 1, 1, 0.5, 0, 0.5, 1, 0.5, 1, 1],
        [1, 0.5, 1, 1, 1, 1, 0.5, 1, 1, 1, 0.5, 0, 1, 0.5, 1, 1],
        [1, 1, 0.5, 1, 0.5, 1, 1, 0.5, 0.5, 1, 1, 1, 0, 1, 0.5, 1],
        [0.5, 1, 1, 0.5, 0.5, 1, 1, 1, 0.5, 1, 0.5, 0.5, 1, 0, 1, 0.5],
        [1, 0.5, 1, 1, 1, 0.5, 0.5, 1, 1, 1, 1, 1, 0.5, 1, 0, 1],
        [1, 1, 1, 1, 0.5, 0.5, 1, 0.5, 1, 1, 1, 1, 1, 0.5, 1, 0]
    ])


# --------------------------------------------DIVERGENCE CALCULATION----------------------------------------

#Calculating the expected genetic divergence between genotypes over time
def calculate_divergence(pedigree_data, initial_state_probs, mutation_rate):
    G = transition_matrix_construction(mutation_rate)
    rate_matrix = logm(G)
    divergence_effects = get_divergence_effects()

    # Storing divergence values calculated for each row of pedigree data
    divergences = []

    for row in pedigree_data:
        time_0 = int(row[0])
        time_1 = int(row[1])
        time_2 = int(row[2])

        prob_generation0 = initial_state_probs.copy()

        # Storing genotype probabilities at time1 and time2 for each 16 genotypes
        prob_generation1 = []
        prob_generation2 = []

        time_diff_1_0 = time_1 - time_0
        time_diff_2_0 = time_2 - time_0

        for i in range(16):
            start = np.eye(16)[i]
            prob_gen1 = np.dot(start, expm(time_diff_1_0 * rate_matrix))
            prob_gen2 = np.dot(start, expm(time_diff_2_0 * rate_matrix))

            prob_generation1.append(prob_gen1)
            prob_generation2.append(prob_gen2)

        prob_generation1 = np.array(prob_generation1)
        prob_generation2 = np.array(prob_generation2)

        div = 0
        for j in range(16):
            joint_prob = np.outer(prob_generation1[j], prob_generation2[j])
            weighted_div = prob_generation0[j] * np.sum(divergence_effects * joint_prob)
            div += weighted_div

        divergences.append(div)

    return np.array(divergences)


# --------------------------------------------CALCULATION FOR PLOT----------------------------------------

# Calculating predicted divergence values and delta times for each pedigree row with optimized gamma and intercept
def calculation_plot_values(pedigree_data, initial_state_probs, best_gamma, best_intercept):
    predicted_divergence = calculate_divergence(pedigree_data[:, :3], initial_state_probs, best_gamma)
    total_predicted_divergence = predicted_divergence + best_intercept
    # ∆t = t1 + t2 - 2 * t0
    delta_t_final = pedigree_data[:, 1] + pedigree_data[:, 2] - 2 * pedigree_data[:, 0]

    plot_values = np.column_stack((pedigree_data[:, :3],  # time0, time1, time2 columns
                                   total_predicted_divergence,  # predicted divergence for each row
                                   delta_t_final))  # time differences for each row

    plot_values = plot_values[plot_values[:, 4].argsort()]  # sorting by delta times

    return plot_values

# -------------------------------------------------LSE----------------------------------------------

# Calculating least squares error between observed and predicted divergence values
def LSE(parameters, pedigree_data, initial_state_probs):
    gamma = parameters[0]
    intercept = parameters[1]

    predicted_divergence = calculate_divergence(pedigree_data, initial_state_probs, gamma) + intercept
    observed_divergence = pedigree_data[:, 3]
    lse_value = np.sum((observed_divergence - predicted_divergence) ** 2)

    return lse_value

# -----------------------------------------------MUTSOMA MAIN------------------------------------------
import numpy as np
from scipy.optimize import minimize
import pprint


def mutSoma(pedigree_path, base_probs, prop_het=0.1, num_starts=50):
    pedigree_data = np.loadtxt(pedigree_path, delimiter='\t', skiprows=1)
    #delta_times = pedigree_data[:, 1] + pedigree_data[:, 2] - 2 * pedigree_data[:, 0]
    initial_state_probs = initial_probabilities(base_probs, prop_het)

    # Normalize pedigree so all time0 values are zero
    pedigree_data[:, 1] -= pedigree_data[:, 0]
    pedigree_data[:, 2] -= pedigree_data[:, 0]
    pedigree_data[:, 0] = 0

    # Storing all optimization results
    opt_results = []

    for i in range(num_starts):
        gamma_start = 10 ** np.random.uniform(np.log10(1e-11), np.log10(1e-3))
        intercept_start = np.random.uniform(0, np.max(pedigree_data[:, 3]))
        initial_params = [gamma_start, intercept_start]

        #Optimization
        result = minimize(
            LSE,
            initial_params,
            args=(pedigree_data, initial_state_probs),
            method="Nelder-Mead",
            options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 5000, "disp": False}
        )

        # Store results
        result_info = {
            "gamma": float(result.x[0]),
            "intercept": float(result.x[1]),
            "success": bool(result.success),
            "LSE": float(result.fun),
            "g.start": float(gamma_start),
            "intercept.start": float(intercept_start),
            "run_id": i + 1
        }

        opt_results.append(result_info)

    valid_results = []
    for res in opt_results:
        if res["success"] and res["gamma"] > 0 and res["intercept"] > 0:
            valid_results.append(res)

    valid_results.sort(key=lambda result: result["LSE"])

    flagged_results = []
    for res in opt_results:
        if res not in valid_results:
            flagged_results.append(res)

    # Summary
    if len(valid_results) > 0:
        best_estimations = valid_results[0]

        for_fit_plot = calculation_plot_values(pedigree_data, initial_state_probs, best_estimations["gamma"],
                                               best_estimations["intercept"])

        print("Run       LSE             Gamma        Intercept              Gamma.start      Intercept.start")
        print("=======================================================================================================")

        for res in valid_results:
            nm_label = "Nelder-Mead" + str(res["run_id"])
            print(nm_label, res["LSE"], res["gamma"], res["intercept"], res["g.start"], res["intercept.start"])

        print("=======================================================================================================")

        # Output dictionary
        output = {
            "estimates": valid_results,
            "estimates_flagged": flagged_results,
            "settings": [
                ["p0aa", base_probs[0]],
                ["p0cc", base_probs[1]],
                ["p0tt", base_probs[2]],
                ["p0gg", base_probs[3]],
                ["Nstarts", num_starts],
                ["prop.het", prop_het],
                ["optim.method", "Nelder-Mead"]
            ],
            "model": "mutSOMA.py",
            "for_fit_plot": for_fit_plot.tolist(),
            "input": pedigree_data
        }

        pprint.pprint(output)
        return output

    else:
        print("No valid optimization results were found")
        return {}


# --------------------------------------------PLOTTING--------------------------------------------
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_predicted_vs_observed(for_fit_plot, used_input):
    observed_delta_t = used_input[:, 1] + used_input[:, 2] - 2 * used_input[:, 0]
    observed_div = used_input[:, 3]

    predicted_delta_t = np.array([row[4] for row in for_fit_plot])
    predicted_div = np.array([row[3] for row in for_fit_plot])

    plt.figure(figsize=(12, 7))

    # Plot observed points
    plt.scatter(
        observed_delta_t, observed_div,
        color='royalblue', s=40, label='Observed',
        alpha=0.6, edgecolor='black', linewidth=0.5
    )

    # Plot predicted line
    plt.plot(
        predicted_delta_t, predicted_div,
        color='#E30B5D', linewidth=2.5, label='Predicted'
    )

    plt.scatter(
        predicted_delta_t, predicted_div,
        color='#5d2cc7', s=40, marker='X', label='Predicted Points',
        alpha=0.8, edgecolor='black', linewidth=0.5
    )

    # Styling
    plt.title("Predicted vs Observed Divergence", fontsize=16, weight='bold')
    plt.xlabel("Delta Time (Δt)", fontsize=14)
    plt.ylabel("Divergence", fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(visible=True, linestyle='--', alpha=0.5)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig("/Users/adaakinci/PycharmProjects/mutSOMA2/plot.png", dpi=300)
    plt.close()


# -----------------------------------------------MAIN--------------------------------------------
# Entry point
if __name__ == "__main__":
    import time
    base_probs = [0.3317790, 0.1685361, 0.3312982, 0.1683867]
    pedigree_path = "/Users/adaakinci/PycharmProjects/mutSOMA2/makeVCFpedigree_output.txt"

    start_time = time.time()  # ⏱ start timing

    output = mutSoma(pedigree_path, base_probs)

    end_time = time.time()  # ⏱ end timing
    print(f"⏱ Runtime: {end_time - start_time:.4f} seconds")

    #output = mutSoma(pedigree_path, base_probs)

    best = output["estimates"][0]
    print("Optimized gamma:", best["gamma"])
    print("Optimized intercept:", best["intercept"])
    initial_state_probs = initial_probabilities(base_probs, prop_het=0.1)
    used_input = output["input"]

    div_obs = used_input[:, 3]
    div_pred = calculate_divergence(used_input[:, :3], initial_state_probs, best["gamma"]) + best["intercept"]

    print("First 5 observed divergence values:", div_obs[:28])
    print("First 5 predicted divergence values:", div_pred[:28])

    predicted = calculate_divergence(used_input, initial_state_probs, best["gamma"])
    print("First 5 predicted (no intercept):", predicted[:5])

    plot_predicted_vs_observed(output["for_fit_plot"], output["input"])
