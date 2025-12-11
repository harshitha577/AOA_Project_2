import random
import time
import matplotlib.pyplot as plt
import numpy as np
from collections import deque

class Dinic:
    def __init__(self, N):
        self.N = N
        self.adj = [[] for _ in range(N)]

    class Edge:
        def __init__(self, to, cap, rev):
            self.to = to
            self.cap = cap
            self.rev = rev

    def add_edge(self, u, v, cap):
        fwd = self.Edge(v, cap, None)
        rev = self.Edge(u, 0, fwd)
        fwd.rev = rev
        self.adj[u].append(fwd)
        self.adj[v].append(rev)

    def bfs(self, s, t, level):
        queue = deque([s])
        level[s] = 0
        while queue:
            v = queue.popleft()
            for e in self.adj[v]:
                if e.cap > 0 and level[e.to] < 0:
                    level[e.to] = level[v] + 1
                    queue.append(e.to)
        return level[t] >= 0

    def dfs(self, v, t, flow, level, it):
        if v == t:
            return flow
        for i in range(it[v], len(self.adj[v])):
            it[v] = i
            e = self.adj[v][i]
            if e.cap > 0 and level[e.to] == level[v] + 1:
                pushed = self.dfs(e.to, t, min(flow, e.cap), level, it)
                if pushed > 0:
                    e.cap -= pushed
                    e.rev.cap += pushed
                    return pushed
        return 0

    def max_flow(self, s, t):
        flow = 0
        INF = 10**18
        while True:
            level = [-1] * self.N
            if not self.bfs(s, t, level):
                return flow
            it = [0] * self.N
            while True:
                pushed = self.dfs(s, t, INF, level, it)
                if pushed <= 0:
                    break
                flow += pushed


def build_and_solve(n_tasks, n_servers):
    demands = [random.randint(2, 6) for _ in range(n_tasks)]
    caps = [random.randint(4, 10) for _ in range(n_servers)]

    compat = [[random.choice([0, 1]) for _ in range(n_servers)]
              for _ in range(n_tasks)]

    SRC = 0
    TASK_START = 1
    SERV_START = TASK_START + n_tasks
    SNK = SERV_START + n_servers
    N = SNK + 1

    dinic = Dinic(N)

    for i in range(n_tasks):
        dinic.add_edge(SRC, TASK_START + i, demands[i])

    for i in range(n_tasks):
        for j in range(n_servers):
            if compat[i][j] == 1:
                dinic.add_edge(TASK_START + i, SERV_START + j, 10**9)

    for j in range(n_servers):
        dinic.add_edge(SERV_START + j, SNK, caps[j])

    total_demand = sum(demands)

    start = time.time()
    maxflow = dinic.max_flow(SRC, SNK)
    end = time.time()

    return maxflow, total_demand, end - start, n_tasks * n_servers


task_sizes = [10, 20, 30, 40, 50, 60]
runtimes = []
edges_list = []

for n in task_sizes:
    mf, td, rt, edges = build_and_solve(n, n)
    runtimes.append(rt)
    edges_list.append(edges)

# Convert for plotting
n_arr = np.array(task_sizes)
rt_arr = np.array(runtimes)
edges_arr = np.array(edges_list)


plt.figure(figsize=(7, 5))
plt.plot(n_arr, rt_arr, marker='o')
plt.xlabel("Number of Tasks / Servers (n = m)")
plt.ylabel("Runtime (seconds)")
plt.title("Runtime Scaling of Max-Flow Reduction")
plt.grid(True)
plt.savefig("runtime_vs_tasks.png", dpi=300)
plt.show()


plt.figure(figsize=(7, 5))
plt.plot(edges_arr, rt_arr, marker='o')
plt.xlabel("Number of Edges in Compatibility Graph")
plt.ylabel("Runtime (seconds)")
plt.title("Runtime vs Number of Edges")
plt.grid(True)
plt.savefig("runtime_vs_edges.png", dpi=300)
plt.show()


plt.figure(figsize=(7, 5))
plt.loglog(n_arr, rt_arr, marker='o')
plt.xlabel("log(Number of Tasks/Servers)")
plt.ylabel("log(Runtime)")
plt.title("Log–Log Runtime Growth of Max-Flow Algorithm")
plt.grid(True, which="both")
plt.savefig("loglog_runtime_plot.png", dpi=300)
plt.show()
