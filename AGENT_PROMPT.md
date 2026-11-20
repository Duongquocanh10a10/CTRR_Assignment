# ROLE
You are a senior C++ engineer and algorithms TA. Implement a university assignment (HCMUT, CO1007 Discrete Structure, "A* Algorithms") in C++, with tests and a LaTeX report. Work in phases and stop for review at the end of each phase.

# HARD RULES (violating any of them costs the student points)
- R1. Language: C++ (assume C++11 or later). The ONLY allowed includes in submitted code: `<iostream> <fstream> <string> <cmath> <vector> <algorithm>`. Do NOT use `<queue> <set> <map> <unordered_map> <limits> <climits> <cstdlib> <cstdio> <utility> <functional> <memory>`. Implement your own binary heap on `std::vector`. Use a large constant (e.g. 1e18) for infinity. Test code under `tests/` may use any standard header, but the submission folder may not.
- R2. Exact function names and signatures (the auto-tester calls them directly):
```cpp
PathNode* findMeshPath(double weightMatrix[100][100], int startNode, int goalNode);
PathNode* findAUVPath(double weightMatrix[100][100], int coords[100][3], int startPoint, int goalPoint, int mode);
PathNode* findProbabilisticPath(double distMatrix[100][100], double riskMatrix[100][100], int coords[100][2], int startPoint, int goalPoint, double lambda, int mode);
PathNode* findGamePath(int mapGrid[100][100], int m, int n, int startX, int startY, int goalX, int goalY, int mode);
PathNode* findPCBPath(int pcbGrid[100][100], int m, int n, int startX, int startY, int exitX, int exitY, double weightMatrix[100][100], int mode);
```
- R3. The common struct must be exactly this (no extra fields):
```cpp
struct PathNode { std::string name; double f; double g; double h; PathNode* next; };
```
- R4. Do NOT define `printPath` in any submitted file (`PathNode.*`, `Algo.*`). The LMS provides it and it prints name/f/g/h. `main.cpp` may contain its own demo copy named `printPath` (main.cpp is not sent to CodeRunner).
- R5. No `using namespace std;` in headers. Use `std::` explicitly in `.h/.cpp`.
- R6. Return value: pointer to the head of a singly linked list ordered start -> goal. No path -> `nullptr`. start == goal -> a single node (g = 0, f = h = 0). Start or goal on an obstacle -> `nullptr`. Each node stores that node's final g, its h, and f = g + h.
- R7. Compare doubles with an epsilon of 1e-9.
- R8. Tie-breaking lives in ONE comparator function shared by all tasks: smaller f first, then smaller h, then smaller node index (or cell index). Isolate it so it can be changed in one line.
- R9. Node names are produced by ONE function per task, `formatName(...)`, so the naming scheme can be changed later without touching the search logic. Defaults: Tasks 1-3 use the vertex index as a string ("0", "1", ...). Task 4 uses direction labels, with the first node named "Start". Task 5 uses "(x,y)".
- R10. The only global variable allowed is `double g_lastPsafe` (Task 3). Declare it `extern` in `Algo.h` and define it in `Algo.cpp`.
- R11. Never invent results. Every number in the report must come from a program run in this environment. If you cannot run something, say so.
- R12. Do not claim a heuristic is admissible unless you proved it. State clearly when it is not.
- R13. Free memory: provide `void freePath(PathNode* head)` in `PathNode.cpp`; use it in tests and `main.cpp`.

# ASSUMPTIONS (the spec is silent; implement them, but keep each one a named constant or a tiny isolated function, and list every one in `NOTES.md`)
- A1. Tasks 1-3 have no node-count parameter: treat the node count as 100. A node with no edges is harmless. Assume the matrices are zero-initialized.
- A2. Edge u->v exists iff `matrix[u][v] > 0` (directed, using the matrix as given). Ignore `u == v`.
- A3. Task 1 `degree(v)` = number of entries `> 0` in ROW v (out-degree). Put it in a single function `degreeOf(...)` so it can be swapped for an in+out or symmetric definition.
- A4. Grid coordinate convention (Tasks 4, 5): `x` is the ROW index (`grid[x][y]`), `y` is the column; Up = x-1, Down = x+1, Left = y-1, Right = y+1. Put this in one place.
- A5. Task 3 modes: 1 = Manhattan, 2 = Euclidean.
- A6. Task 5 grid values: 0 = free, 1 = obstacle. Movement is 4-directional. Node id of cell (x,y) = `x*n + y`. Because the weight matrix is only 100x100, require `m*n <= 100`; otherwise return `nullptr`.
- A7. Allow re-opening: when a cheaper g is found for a node already closed, update it and re-insert it (use lazy deletion in the heap). Stop when the goal is popped from the heap.

# GENERIC ENGINE
Write one reusable A* core (heap on vector, `best[]` g-table, parent array, lazy deletion, comparator from R8) and call it from all five tasks, with a task-specific neighbour/cost/heuristic callback or a small per-task adapter. No function pointers from `<functional>`; plain function pointers or templates are fine. Path reconstruction follows the parent array from goal back to start, then builds the list in start -> goal order.

# TASK SPECIFICATIONS

## Task 1: Mesh network routing with congestion penalty (`findMeshPath`)
- Step cost u->v: `weightMatrix[u][v] + 0.1 * degree(v)` (the penalty applies to the TARGET node; the start node has none). g(start) = 0.
- h(n): 0 if n is the goal; otherwise `minStep` = the minimum of `weightMatrix[u][v] + 0.1*degree(v)` over all edges in the graph (precompute once). This h is admissible and consistent; document why.
- Names: vertex index as string.

## Task 2: AUV navigation in 3D (`findAUVPath`)
- Step cost: `w(u,v) * (1 + |z_new - z_old| / 100)`, where z is `coords[i][2]` (use floating-point division, not integer division).
- h(n) by `mode`: 1 = 3D Manhattan `|dx|+|dy|+|dz|`; 2 = 3D Euclidean; 3 = 3D Chebyshev `max(|dx|,|dy|,|dz|)`. Target = goal station. An invalid mode falls back to mode 1.
- WARNING: these heuristics are not guaranteed admissible against these edge costs. Do not "fix" them; implement exactly as specified, and let the report discuss it.

## Task 3: Probabilistic pathfinding under hazard risk (`findProbabilisticPath`)
- Step cost: `distMatrix[u][v] + lambda * riskMatrix[u][v]`. Edge exists iff `distMatrix[u][v] > 0`. `riskMatrix[u][v]` is in [0, 1).
- h(n): per `mode` (A5), computed from `coords[i][0..1]` to the goal.
- Output: the linked list as usual, plus the global `g_lastPsafe` = product of `(1 - riskMatrix[u][v])` over the returned path (1.0 for a single-node path, 0.0 if there is no path). Reset it at the start of every call.

## Task 4: Game navigation with multi-terrain obstacles (`findGamePath`)
- Cell values: 0 = flat (cost 1.0), 2 = swamp (cost 2.0), 1 = obstacle (cannot be entered).
- 8 directions. Move cost = terrain cost of the TARGET cell; multiply by 1.5 if the move is diagonal.
- Corner-cutting rule: a diagonal move from (x,y) to (x+dx,y+dy) is valid only if BOTH cardinal neighbours (x+dx,y) and (x,y+dy) are inside the grid and are NOT obstacles (value 1). Check bounds before reading the grid.
- h(n): mode 1 = Manhattan, mode 2 = Chebyshev. Note for the report: with diagonal cost 1.5, Manhattan can overestimate (true cost of a 1-cell diagonal is 1.5, Manhattan says 2), so it is not admissible; Chebyshev is admissible.
- Direction labels, exactly: `Up, Down, Left, Right, Up-Left, Up-Right, Down-Left, Down-Right`. The first node is named by `formatName` as "Start"; every later node is named by the direction taken to reach it.

## Task 5: VLSI PCB wire routing with bending penalty (`findPCBPath`)
- First fill the provided `weightMatrix` (zero the used `m*n` x `m*n` block first): weight 1 between 4-adjacent cells that are both free (symmetric), 0 otherwise. This is the "grid -> graph" transformation required by the report.
- Search state = (cell, previous direction) with an extra "none" direction for the start. This is required for correctness: g depends on the previous move. Accept the goal in any direction.
- g: +1 per step; add the bend penalty B = 0.5 when the current move direction differs from the previous direction. The first move has no penalty.
- h(n): mode 1 = Manhattan, mode 2 = Chebyshev (cell to exit); independent of direction.
- Names: "(x,y)" via `formatName`.

# PROJECT LAYOUT
```
submission/            <- ONLY these 5 files; this folder is later renamed to the student ID
  PathNode.h  PathNode.cpp
  Algo.h      Algo.cpp
  main.cpp
tests/                 <- outside the submission folder
  reference/           <- independent brute-force solvers (plain Dijkstra over the same state space, no heuristic)
  test_*.cpp  run_tests.sh
tools/
  make_codrunner.sh    <- for each task, concatenate PathNode + the needed part of Algo into ONE single-file source with no local #include and no printPath, for pasting into CodeRunner (we do not know whether CodeRunner takes multiple files)
report/                <- LaTeX (English), written but not forced to compile if no TeX is installed
NOTES.md               <- assumptions, decisions, open questions
```
`main.cpp` must demonstrate all 5 tasks on small sample inputs (including the Fig. 1 example S->A->D->G with h(S)=3, h(A)=2, h(B)=2, h(C)=1, h(D)=1, h(G)=0 and edges S-A 1, S-B 4, A-C 2, A-D 2, B-D 1, C-G 3, D-G 1, whose optimal cost is 4).

# TESTING
- For every task, compare against the brute-force reference on many random small instances: whenever the heuristic is admissible for that instance, A* cost must equal the reference cost. When it is not, record the gap instead of asserting equality.
- Edge cases for every task: start == goal, no path, start/goal on an obstacle, tie situations, the maximum 100-node size, zero-weight "no edge" entries, a swamp vs detour trade-off (Task 4), a risk vs length trade-off for several values of lambda including 0 (Task 3), and a case where bending makes a longer straight route win (Task 5).
- Check memory with a leak check if available (`valgrind` or `-fsanitize=address,undefined`). Compile with `-Wall -Wextra`, with zero warnings.
- Verify the submission files compile with ONLY the six allowed headers (grep the includes).

# REPORT (LaTeX, English, `report/main.tex`; the final file will be named `<StudentID>.pdf`)
Per task, cover exactly what the assignment asks and use REAL outputs from the code:
1. Mesh: how degrees are computed from the adjacency matrix; how h(n) is formulated and why it is admissible; a manual step-by-step trace on a small network, matching the program output.
2. AUV: compare the three 3D heuristics on concrete test cases: nodes expanded, runtime, path cost, optimality vs the reference. Discuss admissibility.
3. Hazard risk: the effect of lambda (short/risky vs long/safe) with a worked numeric example, including P_safe.
4. Game grid: a step-by-step trace on a small grid showing expansion order, corner-cutting checks and cost updates.
5. PCB: the grid-to-graph algorithm, and how the previous-direction state is kept during the search so bending penalties are computed correctly.
Also add a short section listing the assumptions from NOTES.md.

# STUDY NOTES (important)
The student must be able to explain every line during an oral check, and academic rules forbid unadapted AI output. So:
- Comment the code with the REASON for each non-obvious step, in plain language.
- Create `STUDY_GUIDE.md` that explains, task by task, the algorithm, the key invariants, and 3 likely oral-exam questions with answers.
- Keep the code straightforward and readable rather than clever.

# WORKFLOW
Phase 0: read this prompt, write `NOTES.md` with the assumptions and any questions; do not write code yet. Phase 1: `PathNode` + engine + Task 1 + tests. Phase 2: Tasks 2 and 3. Phase 3: Tasks 4 and 5. Phase 4: `main.cpp`, `make_codrunner.sh`, compile-check on the allowed headers only. Phase 5: report and STUDY_GUIDE. At the end of each phase, stop and print: what was done, test results (real), and open questions. If you find a conflict inside this prompt or a spec ambiguity not covered above, do not silently guess: add it to NOTES.md under "Open questions" and pick the most conservative option.

# DEFINITION OF DONE
All five functions compile with only the six headers and the exact signatures; all tests pass or have recorded, explained gaps; `main.cpp` runs the demos; NOTES.md lists every assumption; the report contains only real numbers.
