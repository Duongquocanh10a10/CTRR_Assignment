# STUDY GUIDE — CO1007 Discrete Structure Assignment (A* Algorithms)

This guide summarizes the shared A* architecture, each task's algorithm and invariants, and likely oral-exam questions with clear answers so you can explain every line of code confidently.

---

## 0. Shared Engine (`PathNode.*` and top of `Algo.cpp`)

- **Data structure:** `PathNode { std::string name; double f; double g; double h; PathNode* next; }`.
- **Why a custom binary min-heap?** LMS CodeRunner only allows `<iostream>`, `<fstream>`, `<string>`, `<cmath>`, `<vector>`, `<algorithm>` — `<queue>` is forbidden. `BinaryMinHeap` stores `HeapEntry {f, g, h, state}` in a `std::vector` with `siftUp` / `siftDown` in $O(\log N)$ time.
- **Why lazy deletion?** Binary heaps don't support $O(\log N)$ `decrease-key` without an index map. Instead, whenever we find a cheaper `g` (`tentativeG < bestG[nextState] - EPS`), we push a new entry and later skip any popped entry whose `cur.g > bestG[cur.state] + EPS`.
- **Tie-breaking comparator (`isHigherPriority`):** Smaller `f` (beyond `EPS = 1e-9`), then smaller `h`, then smaller `state` index. Isolated in ONE function so it can be changed in 1 line if LMS test cases use a different priority.

---

## 1. Task 1 — Mesh Network Packet Routing (`findMeshPath`)

- **Algorithm:**
  1. Precompute `degrees[v] = degreeOf(weightMatrix, v)` (out-degree: count of `j != v` with `weightMatrix[v][j] > 0`).
  2. Precompute `minStep` = minimum of `weightMatrix[u][v] + 0.1 * degrees[v]` over all edges `u -> v` (`0.0` if no edges).
  3. Run `runAStar` with step cost `weightMatrix[u][v] + 0.1 * degrees[v]` and heuristic `h(n) = (n == goal ? 0.0 : minStep)`.
- **Key invariants:**
  - Congestion penalty applies to the **target** node `v`, never the start node (`g(start) = 0`).
  - `h(n)` is both **admissible** (any non-empty path takes at least 1 edge of cost $\ge \text{minStep}$) and **consistent** ($h(u) \le c(u,v) + h(v)$ for every edge).
- **Likely oral-exam questions:**
  1. *Q: Why is `h(n) = minStep` admissible when we have no coordinates?*
     **A:** If $n \neq \text{goal}$, any path from $n$ to $\text{goal}$ must traverse at least one edge $(u,v)$, and every edge costs at least `minStep`. Thus $h(n) = \text{minStep} \le h^*(n)$. At the goal, $h(\text{goal}) = 0 = h^*(\text{goal})$.
  2. *Q: Why did you precompute `degrees[]` instead of recomputing inside `forEachNeighbor`?*
     **A:** Precomputing takes $O(V^2)$ once. Recomputing inside the neighbor loop would add an extra $O(V)$ scan per relaxed edge, making the search $O(V^3)$ instead of $O((V + E)\log V)$.
  3. *Q: What happens if a graph is undirected vs. directed for `degree(v)`?*
     **A:** On an undirected (symmetric) adjacency matrix, row count (out-degree) equals column count (in-degree) equals the undirected degree. Isolating `degreeOf` in one function lets us swap the definition cleanly if needed.

---

## 2. Task 2 — AUV Submarine Navigation in 3D (`findAUVPath`)

- **Algorithm:**
  - Step cost `u -> v`: `weightMatrix[u][v] * (1.0 + fabs(z_v - z_u) / 100.0)` using `double` division (not integer division!).
  - Heuristic `h(n)` to `goalPoint`: Mode 1 = 3D Manhattan, Mode 2 = 3D Euclidean, Mode 3 = 3D Chebyshev (invalid mode falls back to Mode 1).
- **Key invariants:**
  - Always $h_{\text{Chebyshev}} \le h_{\text{Euclidean}} \le h_{\text{Manhattan}}$.
  - **Not unconditionally admissible:** Edge weights `w(u,v)` are arbitrary energy costs, not guaranteed to be $\ge$ spatial distance.
- **Likely oral-exam questions:**
  1. *Q: Why must we cast `coords[v][2] - coords[u][2]` to `double` before dividing by `100.0`?*
     **A:** `coords` is `int[100][3]`. Dividing an `int` by `100` in C++ performs integer truncation (e.g. `99 / 100 == 0`), silently wiping out depth penalties smaller than 100 units. Dividing by `100.0` (with `static_cast<double>`) preserves exact fractional penalties.
  2. *Q: Are the 3D heuristics admissible here?*
     **A:** Not in general. If an edge has low energy cost `w(u,v) = 1` across a 50-unit coordinate jump, $h(u) = 50 > c(u,v) = 1$, overestimating true cost. They ARE admissible whenever $w(u,v) \ge d_k(u,v)$ for every edge, since $(1 + |\Delta z|/100) \ge 1$.
  3. *Q: When all three heuristics are admissible, which one expands the fewest nodes and why?*
     **A:** 3D Manhattan expands the fewest nodes because it dominates the others ($h_{\text{Manhattan}} \ge h_{\text{Euclidean}} \ge h_{\text{Chebyshev}}$), giving tighter $f$-values that prune more off-path nodes.

---

## 3. Task 3 — Probabilistic Pathfinding under Hazard Risk (`findProbabilisticPath`)

- **Algorithm:**
  - Reset `g_lastPsafe = 0.0` immediately at function entry.
  - Step cost `u -> v`: `distMatrix[u][v] + lambda * riskMatrix[u][v]` (edge exists iff `distMatrix[u][v] > 0`).
  - Heuristic `h(n)`: Mode 1 = 2D Manhattan, Mode 2 = 2D Euclidean on `coords[i][0..1]`.
  - After `runAStar` returns `pathStates`, compute `g_lastPsafe = prod(1.0 - riskMatrix[u][v])` (`1.0` for single-node path, `0.0` if no path).
- **Key invariants:**
  - `g` is additive (expected cost), whereas `P_safe` is multiplicative (assuming independent edge survival events).
  - Increasing $\lambda$ shifts preference from short/risky routes toward longer/safer routes.
- **Likely oral-exam questions:**
  1. *Q: Why don't we maximize `P_safe` directly inside A* instead of using `lambda * riskMatrix[u][v]`?*
     **A:** The assignment defines the objective as a linear expected-cost trade-off $g(n) = \sum [d(u,v) + \lambda P_{\text{risk}}(u,v)]$ so $\lambda$ controls the balance between distance and risk, and then evaluates $P_{\text{safe}}$ on the chosen route.
  2. *Q: Why do we reset `g_lastPsafe = 0.0` at the start of every call?*
     **A:** Because `g_lastPsafe` is a global variable, a previous successful call could leave a stale non-zero value if the next call fails (unreachable or invalid indices) unless we reset it upfront.
  3. *Q: What if `startPoint == goalPoint`?*
     **A:** The path has zero edges, so the empty product gives $P_{\text{safe}} = 1.0$, and the returned linked list has a single node with $g = h = f = 0$.

---

## 4. Task 4 — Game Character Navigation (`findGamePath`)

- **Algorithm:**
  - Grid of size `m x n` (`x` = row, `y` = col), 8 directions (`DIR_DX`, `DIR_DY`, `DIR_NAMES`).
  - Terrain cost of **target** cell: `0` (flat) = `1.0`, `2` (swamp) = `2.0`, `1` (obstacle) = blocked.
  - Diagonal move multiplies base cost by `1.5`, and is allowed ONLY if both cardinal neighbors `(x + dx, y)` and `(x, y + dy)` are in bounds and NOT `1` (corner-cutting rule).
  - Heuristic: Mode 1 = Manhattan, Mode 2 = Chebyshev.
  - Node naming: `"Start"` for first node, then direction taken to reach each subsequent node.
- **Key invariants:**
  - Mode 2 (Chebyshev) is admissible and consistent (every step costs $\ge 1.0$ and changes Chebyshev distance by at most 1).
  - Mode 1 (Manhattan) is **not admissible** because a 1-step flat diagonal costs `1.5` while Manhattan estimates `2.0`.
- **Likely oral-exam questions:**
  1. *Q: Why is Manhattan not admissible in Task 4, whereas Chebyshev is?*
     **A:** Moving diagonally from `(0,0)` to `(1,1)` on flat ground has true cost `1.5`, but Manhattan gives $|1-0| + |1-0| = 2.0 > 1.5$ (overestimates). Chebyshev gives $\max(1, 1) = 1.0 \le 1.5$, and no move can ever reduce $\max(|\Delta x|, |\Delta y|)$ by more than 1 while costing less than `1.0`.
  2. *Q: How does your code prevent diagonal corner-cutting?*
     **A:** Before accepting a diagonal move `(x, y) -> (x + dx, y + dy)`, we check bounds and verify `mapGrid[x + dx][y] != 1 && mapGrid[x][y + dy] != 1`. Note that swamp cells (`2`) do NOT block corner-cutting — only obstacles (`1`) do.
  3. *Q: How do you reconstruct the direction labels from the parent array?*
     **A:** `runAStar` passes `(curState, prevState)` to `adapter.formatName`. If `prevState == -1`, it returns `"Start"`; otherwise it computes `(dx, dy)` between `prevCell` and `curCell` and looks up the matching string in `DIR_NAMES`.

---

## 5. Task 5 — VLSI PCB Wire Routing with Bending Penalty (`findPCBPath`)

- **Algorithm:**
  1. Check `m > 0`, `n > 0`, `m * n <= 100`. Zero the `(m*n) x (m*n)` block of `weightMatrix`, then set `weightMatrix[u][v] = 1.0` for every 4-adjacent pair of free cells (`pcbGrid == 0`).
  2. Validate `start` and `exit` are in bounds and free.
  3. Run `runAStar` on the augmented state space `state = cellId * 5 + prevDir` (`prevDir` in `0..3` = Up/Down/Left/Right, `4` = None for start).
  4. Step cost = `1.0 + (prevDir != 4 && prevDir != d ? 0.5 : 0.0)`. Goal test = `(state / 5) == exitCell`.
  5. Heuristic: Mode 1 = Manhattan, Mode 2 = Chebyshev (independent of direction; both admissible and consistent).
- **Key invariants:**
  - Augmenting the state with `prevDir` is **required for correctness** — without it, a cheaper arrival to an intermediate cell in the "wrong" direction could prune a slightly more expensive arrival in the "right" direction that avoids a future `0.5` bend.
  - Because `weightMatrix` is `100 x 100`, `m * n` must be $\le 100$.
- **Likely oral-exam questions:**
  1. *Q: Why can't we just store `parent[cell]` and check direction on the fly without expanding the state space to `m * n * 5`?*
     **A:** Suppose two paths reach cell `C` with the same cost `g = 4.0`: Path 1 arrives going Right, Path 2 arrives going Down. If the only exit from `C` goes Down, Path 2 continues for `1.0` while Path 1 pays `1.5`. If we only kept one `bestG[C]`, Path 1 would discard Path 2 and return a suboptimal path! Separating `(C, Right)` and `(C, Down)` keeps both.
  2. *Q: Why is `prevDir = 4` (`None`) needed for the start state?*
     **A:** The first wire segment leaving the start pin has no prior direction, so it must NOT incur a `0.5` bending penalty regardless of whether it starts Up, Down, Left, or Right.
  3. *Q: Why are both Manhattan and Chebyshev admissible in Task 5?*
     **A:** Movement is 4-directional and every edge costs at least `1.0` (possibly `1.5` with a bend). Since each step reduces Manhattan distance by at most `1`, we have $\text{Chebyshev} \le \text{Manhattan} \le \text{true cost}$.
