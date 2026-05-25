import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

np.random.seed(42)

# =====================================================
# PARAMETERS
# =====================================================

GENERATIONS = 250

PEAK_POP = 20000
PEAK_TIME = 120
EPIDEMIC_WIDTH = 30 #45

MUTATION_RATE = 0.015

INITIAL_T = 1.0
INITIAL_E = 0.02

TRADEOFF = 0.25

BASE_MUT_STD = 0.03
ESCAPE_SCALE = 0.08

SHIFT_POINT = 175
K = 0.07

# =====================================================
# EPIDEMIC CURVE  ¡¾change here for constant population¡¿
# =====================================================

def population_size(t):

    # bell-shaped epidemic curve
    size = PEAK_POP * np.exp(
        -((t - PEAK_TIME) ** 2) /
        (2 * EPIDEMIC_WIDTH ** 2)
    )

    return max(int(size), 200)
    #if t <= PEAK_TIME:
    #    return max(int(size), 200)
    #else:
    #    return 20000

# =====================================================
# SELECTION SHIFT
# =====================================================

def immune_weight(t):

    return 1 / (1 + np.exp(-K * (t - SHIFT_POINT)))

# =====================================================
# HAPLOTYPE
# =====================================================

class Haplotype:

    def __init__(self, T, E, lineage):

        self.T = T
        self.E = E
        self.lineage = lineage

# =====================================================
# FITNESS  ¡¾change here for constant selection¡¿
# =====================================================

def fitness(h, t):

    beta = immune_weight(t)
    alpha = 1 - beta

    #beta = 0.5
    #alpha = 0.5

    return alpha * h.T + beta * h.E

# =====================================================
# INITIAL POPULATION
# =====================================================

population = []

for i in range(500):

    population.append(
        Haplotype(
            INITIAL_T,
            INITIAL_E,
            "root"
        )
    )

# =====================================================
# TRACKING
# =====================================================

mean_E = []
mean_T = []
max_escape = []
lag_E = []
diversity = []
pop_sizes = []
dominant_E = []
alpha_history = []
beta_history = []

lineage_counter = 0

# =====================================================
# SIMULATION
# =====================================================

for t in range(GENERATIONS):

    target_pop = population_size(t)

    fit = np.array([
        fitness(h, t)
        for h in population
    ])

    beta_history.append(immune_weight(t))
    alpha_history.append(1-immune_weight(t))
    #beta_history.append(0.5)
    #alpha_history.append(0.5)

    fit = fit - fit.min() + 0.001

    offspring_probs = fit / fit.sum()

    # select parents
    parent_indices = np.random.choice(
        np.arange(len(population)),
        size=target_pop,
        p=offspring_probs
    )

    new_population = []

    for idx in parent_indices:

        parent = population[idx]

        child_T = parent.T
        child_E = parent.E
        lineage = parent.lineage

        # mutation
        if np.random.rand() < MUTATION_RATE:

            lineage_counter += 1
            lineage = f"L{lineage_counter}"

            # ------------------------------------------------
            # Most mutations:
            # small/no escape gain
            # occasional larger jumps
            # ------------------------------------------------

            delta_E = np.random.gamma(
                shape=1.2,
                scale=ESCAPE_SCALE
            )
            #delta_E = np.random.normal(
            #loc=0,
            #scale=0.03
            #)

            # many mutations have tiny effect
            delta_E *= np.random.binomial(1, 0.25)

            # transmission tradeoff
            delta_T = (
                np.random.normal(0, BASE_MUT_STD)
                - TRADEOFF * delta_E
            )

            child_E += delta_E

            # prevent negative escape
            child_E = max(child_E, 0)

            child_T += delta_T

        new_population.append(
            Haplotype(
                child_T,
                child_E,
                lineage
            )
        )

    population = new_population

    # =================================================
    # RECORD STATS
    # =================================================

    E_vals = [h.E for h in population]
    T_vals = [h.T for h in population]

    mean_E.append(np.mean(E_vals))
    mean_T.append(np.mean(T_vals))
    max_escape.append(np.max(E_vals))
    lag_E.append(np.max(E_vals)-np.mean(E_vals))
    lineage_counts = defaultdict(int)

    for h in population:
        lineage_counts[h.lineage] += 1

    diversity.append(len(lineage_counts))
    pop_sizes.append(len(population))

    dominant_lineage = max(
        lineage_counts,
        key=lineage_counts.get
    )

    dominant_haps = [
        h for h in population
        if h.lineage == dominant_lineage
    ]

    dominant_E.append(
        np.mean([h.E for h in dominant_haps])
    )

# =====================================================
# PLOTS
# =====================================================

fig, axes = plt.subplots(4, 1, figsize=(9, 12))

axes[0].plot(pop_sizes)
axes[0].axvline(SHIFT_POINT, linestyle='--')
axes[0].set_ylabel("Population size")

axes[1].plot(alpha_history, label="Transmission weight")
axes[1].plot(beta_history, label="Immune escape weight")
axes[1].legend()
axes[1].set_ylabel("Selection weights")

axes[2].plot(dominant_E)
axes[2].plot(max_escape)
axes[2].axvline(SHIFT_POINT, linestyle='--')
axes[2].legend(["Dominant escape","Max escape potential"])
axes[2].set_ylabel("Immune escape")

axes[3].plot(lag_E)
axes[3].axvline(SHIFT_POINT, linestyle='--')
axes[3].legend(["Escape lag"])
axes[3].set_ylabel("Immune escape")
axes[3].set_xlabel("Generation")

plt.tight_layout()
plt.show()
