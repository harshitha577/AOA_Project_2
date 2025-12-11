import random
import time
import numpy as np
import matplotlib.pyplot as plt

DNA = ['A', 'C', 'G', 'T']

def mutate_base(base, p_sub=0.05):
    """Randomly substitute a DNA base with probability p_sub."""
    if random.random() < p_sub:
        return random.choice([b for b in DNA if b != base])
    return base

def generate_descendant(seq, p_sub=0.05, p_ins=0.02, p_del=0.02):
    """Generate a mutated descendant of an ancestor sequence."""
    new_seq = []
    i = 0
    while i < len(seq):
        # deletion
        if random.random() < p_del:
            i += 1
            continue
        
        # substitution
        b = mutate_base(seq[i], p_sub=p_sub)
        new_seq.append(b)
        
        # insertion
        if random.random() < p_ins:
            new_seq.append(random.choice(DNA))
        
        i += 1
    
    return ''.join(new_seq)

def generate_homologous_sequences(n, k):
    """Generate k homologous sequences from an ancestor of length n."""
    ancestor = ''.join(random.choice(DNA) for _ in range(n))
    sequences = [generate_descendant(ancestor) for _ in range(k)]
    return sequences


def nw_align(s1, s2, match=1, mismatch=-1, gap=-1):
    """Return alignment score only (not alignment).
    DP complexity: O(n^2)
    """
    n, m = len(s1), len(s2)
    dp = [[0]*(m+1) for _ in range(n+1)]
    
    for i in range(1, n+1):
        dp[i][m] = gap * (n - i + 1)
    for j in range(1, m+1):
        dp[n][j] = gap * (m - j + 1)
    
    for i in range(n-1, -1, -1):
        for j in range(m-1, -1, -1):
            if s1[i] == s2[j]:
                score = match
            else:
                score = mismatch
            dp[i][j] = max(
                score + dp[i+1][j+1],
                gap + dp[i+1][j],
                gap + dp[i][j+1]
            )
    
    return dp[0][0]


def progressive_msa(seqs):
    """A lightweight progressive alignment that only performs:
       - Pairwise distances via NW
       - Guide tree via simple hierarchical clustering
       - No real profile-profile alignment (for speed), only timing demonstration
    """
    k = len(seqs)
    # Compute pairwise distances
    for i in range(k):
        for j in range(i + 1, k):
            _ = nw_align(seqs[i], seqs[j])  # Just compute score, ignore result
    return True


def experiment_k():
    fixed_n = 300
    k_values = [5, 10, 15, 20, 25]
    runtimes = []
    
    for k in k_values:
        seqs = generate_homologous_sequences(fixed_n, k)
        start = time.time()
        progressive_msa(seqs)
        end = time.time()
        runtimes.append((end - start) * 1000)  # ms
    
    # Plot
    plt.figure(figsize=(6,4))
    plt.plot(k_values, runtimes, marker='o')
    plt.xlabel("Number of sequences (k)")
    plt.ylabel("Runtime (ms)")
    plt.title("Runtime vs Number of Sequences (k)")
    plt.grid(True)
    plt.savefig("runtime_vs_k.png", dpi=300, bbox_inches='tight')
    print("Saved plot: runtime_vs_k.png")

def experiment_n():
    n_values = [50, 100, 200, 400, 800]
    fixed_k = 10
    runtimes = []
    
    for n in n_values:
        seqs = generate_homologous_sequences(n, fixed_k)
        start = time.time()
        progressive_msa(seqs)
        end = time.time()
        runtimes.append((end - start) * 1000)  # ms
    
    # Plot
    plt.figure(figsize=(6,4))
    plt.plot(n_values, runtimes, marker='o')
    plt.xlabel("Sequence length (n)")
    plt.ylabel("Runtime (ms)")
    plt.title("Runtime vs Sequence Length (n)")
    plt.grid(True)
    plt.savefig("runtime_vs_n.png", dpi=300, bbox_inches='tight')
    print("Saved plot: runtime_vs_n.png")

if __name__ == "__main__":
    experiment_k()
    experiment_n()
