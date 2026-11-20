# NOTES — CO1007 Discrete Structure Assignment (A* Algorithms)

This document records all hard constraints, design assumptions, and open questions identified during Phase 0 after reading `AGENT_PROMPT.md` and `DS__CO1007__SEM_261___Assignment___A_ (1).pdf`.

---

## 1. Hard Constraints Summary

- **Language & Headers (R1):** C++11 or later. Submitted files (`submission/PathNode.*`, `submission/Algo.*`, `submission/main.cpp`) may ONLY include:
  `<iostream>`, `<fstream>`, `<string>`, `<cmath>`, `<vector>`, `<algorithm>`.
  No `<queue>`, `<set>`, `<map>`, `<unordered_map>`, `<limits>`, `<climits>`, `<cstdlib>`, `<cstdio>`, `<utility>`, `<functional>`, `<memory>`.
  Custom binary min-heap backed by `std::vector`; use `1e18` for infinity.
- **Exact Signatures (R2):**
  - `PathNode* findMeshPath(double weightMatrix[100][100], int startNode, int goalNode);`
  - `PathNode* findAUVPath(double weightMatrix[100][100], int coords[100][3], int startPoint, int goalPoint, int mode);`
  - `PathNode* findProbabilisticPath(double distMatrix[100][100], double riskMatrix[100][100], int coords[100][2], int startPoint, int goalPoint, double lambda, int mode);`
  - `PathNode* findGamePath(int mapGrid[100][100], int m, int n, int startX, int startY, int goalX, int goalY, int mode);`
  - `PathNode* findPCBPath(int pcbGrid[100][100], int m, int n, int startX, int startY, int exitX, int exitY, double weightMatrix[100][100], int mode);`
- **Common Struct (R3):**
  `struct PathNode { std::string name; double f; double g; double h; PathNode* next; };` (exact fields, no extras).
- **No `printPath` in submitted library files (R4):** LMS provides `printPath`. Only `main.cpp` (demo) may define a local `printPath`.
- **Namespace hygiene (R5):** No `using namespace std;` in headers; use `std::` explicitly in `.h` and `.cpp`.
- **Return & Edge-Case Semantics (R6):**
  - Singly linked list ordered `start -> goal`.
  - No path -> `nullptr`.
  - `start == goal` (and valid/non-obstacle) -> single node (`g = 0`, `h = 0`, `f = 0`).
  - Start or goal on an obstacle (or out of bounds) -> `nullptr`.
  - Every node stores its final `g`, its `h`, and `f = g + h`.
- **Floating-Point Comparison (R7):** Epsilon `EPS = 1e-9`.
- **Single Shared Tie-Breaking Comparator (R8):** Smaller `f`, then smaller `h`, then smaller state/node index (isolated in one function for easy 1-line adjustment).
- **Isolated Naming Helpers (R9):** One `formatName(...)` per task:
  - Tasks 1–3: vertex index as string (`"0"`, `"1"`, ...).
  - Task 4: `"Start"` for the first node, then exact direction labels (`Up`, `Down`, `Left`, `Right`, `Up-Left`, `Up-Right`, `Down-Left`, `Down-Right`).
  - Task 5: `"(x,y)"`.
- **Globals (R10):** Only `extern double g_lastPsafe;` (declared in `Algo.h`, defined in `Algo.cpp`). Reset at the start of every `findProbabilisticPath` call (`1.0` for single-node path, `0.0` if no path).
- **Memory Management (R13):** Provide `void freePath(PathNode* head)` in `PathNode.h`/`PathNode.cpp`.
- **Honesty & Admissibility (R11, R12):** All report numbers must come from real runs; never claim admissibility without proof, and explicitly note when a heuristic is non-admissible.

---

## 2. Assumptions (Isolated in Code for Easy Adjustment)

| ID | Scope | Assumption | Isolation Strategy |
|----|-------|------------|--------------------|
| A1 | Tasks 1–3 | Node count is fixed at `MAX_NODES = 100`; matrices are zero-initialized; isolated nodes with no edges are harmless. | Named constant `MAX_NODES = 100`. |
| A2 | Tasks 1–3 | Directed edge `u -> v` exists iff `matrix[u][v] > 0.0` (using the matrix as given). Self-loops (`u == v`) are ignored. | Helper `hasEdge(w)`. |
| A3 | Task 1 | `degree(v)` = out-degree = count of `j != v` in row `v` where `weightMatrix[v][j] > 0.0`. | Single function `degreeOf(weightMatrix, v)`. |
| A4 | Tasks 4–5 | Grid coordinates: `x` is row (`grid[x][y]`), `y` is column; `Up = x-1`, `Down = x+1`, `Left = y-1`, `Right = y+1`. | Single direction table (`DX[]`, `DY[]`, `DIR_NAMES[]`). |
| A5 | Task 3 | Heuristic modes: `1` = 2D Manhattan, `2` = 2D Euclidean; invalid mode falls back to Mode 1. | Single heuristic helper `heuristicProbabilistic(...)`. |
| A6 | Task 5 | `pcbGrid` values: `0` = free, `1` (non-zero) = obstacle. 4-directional movement. Cell `(x, y)` maps to node ID `x * n + y`. Requires `m > 0`, `n > 0`, and `m * n <= 100`; otherwise return `nullptr`. | Helper `cellToId(x, y, n)` and `isFreePCB(...)`. |
| A7 | All tasks | Re-opening allowed: if a cheaper `g` (strictly less by `EPS`) is found for a state, update `best[state]` and push a new entry into the min-heap (lazy deletion on pop). Stop when the goal is popped from the heap. | Shared generic A* engine in `Algo.cpp`. |

---

## 3. Open Questions & Conservative Defaults

These points are not explicitly specified in the PDF (`DS__CO1007__SEM_261___Assignment___A_ (1).pdf`) or have potential LMS-specific nuances. Each has a conservative default implemented so it can be updated in one place if the LMS FAQ or auto-tester clarifies it:

1. **Task 1 — `degree(v)` on directed vs. undirected graphs:**
   The PDF defines `degree(v)` as "number of active connections of target node `v`" without specifying out-degree vs. in-degree vs. undirected union.
   - *Conservative default (A3):* Out-degree (row count of `weightMatrix[v][j] > 0`, `j != v`), isolated in `degreeOf(...)`. On symmetric undirected matrices this equals the undirected degree.

2. **Task 1 — `h(n)` formula:**
   The PDF only states "$h(n)$: Estimated minimum routing latency remaining to reach the destination node."
   - *Conservative default (per `AGENT_PROMPT.md`):* `h(goal) = 0`, and for `n != goal`, `h(n) = minStep = min_{(u,v) \in E} (weightMatrix[u][v] + 0.1 * degree(v))` (or `0.0` if no edges exist). This is provably admissible and consistent because any non-empty path to the goal takes at least one edge of cost `>= minStep`.

3. **Task 3 — Reporting `P_safe`:**
   The PDF lists accumulated `P_safe` as part of the output, while `PathNode` has no `p_safe` field and `findProbabilisticPath` returns `PathNode*`.
   - *Conservative default (R10):* Store it in `double g_lastPsafe` (declared `extern` in `Algo.h`, defined in `Algo.cpp`), reset at the start of every call (`1.0` if `start == goal`, `0.0` if no path).

4. **Tasks 2–5 — Invalid `mode` parameter:**
   The PDF does not specify behavior for out-of-range `mode` values.
   - *Conservative default:* Fall back to Mode 1 (Manhattan) across Tasks 2, 3, 4, and 5.

5. **Task 4 — Start node label & invalid/obstacle coordinates:**
   The PDF does not specify the label for the start node in the direction-label list.
   - *Conservative default (R9/R6):* First node is named `"Start"` via `formatName`; out-of-bounds or obstacle start/goal returns `nullptr`; `start == goal` (non-obstacle) returns a single node `("Start", f=0, g=0, h=0)`.

6. **Task 5 — Populating `weightMatrix` when start/exit is an obstacle:**
   When `m > 0`, `n > 0`, and `m * n <= 100`, if `start` or `exit` is on an obstacle, should `weightMatrix` still be zeroed and populated before returning `nullptr`?
   - *Conservative default:* Yes — always perform the grid-to-graph population of the `m*n x m*n` block first (whenever `m, n` are valid and `m * n <= 100`), then validate start/exit and run A*.

7. **Tasks 1–3 — Out-of-bounds vertex indices:**
   What if `startNode`/`goalNode` or `startPoint`/`goalPoint` is outside `[0, 99]`?
   - *Conservative default:* Return `nullptr` (and set `g_lastPsafe = 0.0` in Task 3).

8. **Tie-breaking priority for LMS CodeRunner:**
   PDF §4.2 notes: *"Any cases that require considering priority, refer to these test cases for the answers."*
   - *Conservative default (R8):* Smaller `f`, then smaller `h`, then smaller node/state index, isolated in one comparator function so it can be tuned in one line once the LMS test cases open.

9. **Student ID & Packaging:**
   Final submission requires a folder `2410053/` containing `PathNode.h`, `PathNode.cpp`, `Algo.h`, `Algo.cpp`, `main.cpp`, and `2410053.pdf`, compressed as `2410053.rar`.
   - *Confirmed by user:* Student ID = `2410053`. All conservative defaults above approved.
