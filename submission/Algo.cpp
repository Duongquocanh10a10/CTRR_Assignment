#include "Algo.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// ============================================================================
// SHARED ENGINE (Constants, Min-Heap, Tie-Breaking Comparator, Generic A*)
// ============================================================================
namespace {

const double INF = 1e18;
const double EPS = 1e-9;
const int MAX_NODES = 100; // A1: Tasks 1-3 use fixed 100x100 matrices.

// A2: An edge u -> v exists iff its weight is strictly positive.
inline bool hasEdge(double weight) {
    return weight > 0.0;
}

struct HeapEntry {
    double f;
    double g;
    double h;
    int state;
};

// R8: Single shared tie-breaking comparator for all tasks.
// Returns true if entry 'a' has higher priority (should be popped before 'b').
// Order: smaller f, then smaller h, then smaller state index.
bool isHigherPriority(const HeapEntry& a, const HeapEntry& b) {
    if (std::fabs(a.f - b.f) > EPS) {
        return a.f < b.f;
    }
    if (std::fabs(a.h - b.h) > EPS) {
        return a.h < b.h;
    }
    return a.state < b.state;
}

// Custom binary min-heap backed by std::vector (R1: <queue> is forbidden).
class BinaryMinHeap {
public:
    bool empty() const {
        return data_.empty();
    }

    void push(const HeapEntry& entry) {
        data_.push_back(entry);
        siftUp(data_.size() - 1);
    }

    HeapEntry pop() {
        HeapEntry top = data_[0];
        data_[0] = data_.back();
        data_.pop_back();
        if (!data_.empty()) {
            siftDown(0);
        }
        return top;
    }

private:
    std::vector<HeapEntry> data_;

    void siftUp(std::size_t idx) {
        while (idx > 0) {
            std::size_t parent = (idx - 1) / 2;
            if (isHigherPriority(data_[idx], data_[parent])) {
                std::swap(data_[idx], data_[parent]);
                idx = parent;
            } else {
                break;
            }
        }
    }

    void siftDown(std::size_t idx) {
        std::size_t n = data_.size();
        while (true) {
            std::size_t left = 2 * idx + 1;
            std::size_t right = 2 * idx + 2;
            std::size_t best = idx;

            if (left < n && isHigherPriority(data_[left], data_[best])) {
                best = left;
            }
            if (right < n && isHigherPriority(data_[right], data_[best])) {
                best = right;
            }
            if (best != idx) {
                std::swap(data_[idx], data_[best]);
                idx = best;
            } else {
                break;
            }
        }
    }
};

// Generic A* search engine shared across all tasks.
// Supports re-opening via lazy deletion (A7) and stops when a goal state is popped.
template <typename Adapter>
PathNode* runAStar(int numStates,
                   int startState,
                   const Adapter& adapter,
                   std::vector<int>* outStates = nullptr) {
    if (numStates <= 0 || startState < 0 || startState >= numStates) {
        return nullptr;
    }

    std::vector<double> bestG(numStates, INF);
    std::vector<double> bestH(numStates, 0.0);
    std::vector<int> parent(numStates, -1);

    double startH = adapter.heuristic(startState);
    bestG[startState] = 0.0;
    bestH[startState] = startH;

    BinaryMinHeap openSet;
    openSet.push({startH, 0.0, startH, startState});

    int goalState = -1;

    while (!openSet.empty()) {
        HeapEntry cur = openSet.pop();

        // Lazy deletion: ignore stale heap entries that have a worse g than bestG.
        if (cur.g > bestG[cur.state] + EPS) {
            continue;
        }

        // First time a goal state is popped, we have our path (A7).
        if (adapter.isGoal(cur.state)) {
            goalState = cur.state;
            break;
        }

        // Expand outgoing transitions from cur.state.
        adapter.forEachNeighbor(cur.state, [&](int nextState, double stepCost) {
            double tentativeG = cur.g + stepCost;
            // Relax if strictly cheaper (by EPS) than any previously known path.
            if (tentativeG < bestG[nextState] - EPS) {
                double nextH = adapter.heuristic(nextState);
                bestG[nextState] = tentativeG;
                bestH[nextState] = nextH;
                parent[nextState] = cur.state;
                openSet.push({tentativeG + nextH, tentativeG, nextH, nextState});
            }
        });
    }

    if (goalState == -1) {
        return nullptr;
    }

    // Reconstruct state sequence from goal back to start, then reverse to start -> goal.
    std::vector<int> pathStates;
    for (int s = goalState; s != -1; s = parent[s]) {
        pathStates.push_back(s);
    }
    std::reverse(pathStates.begin(), pathStates.end());

    if (outStates != nullptr) {
        *outStates = pathStates;
    }

    // Materialize the PathNode singly linked list in start -> goal order (R6).
    PathNode* head = nullptr;
    PathNode* tail = nullptr;
    for (std::size_t i = 0; i < pathStates.size(); ++i) {
        int curState = pathStates[i];
        int prevState = (i == 0) ? -1 : pathStates[i - 1];
        double gVal = bestG[curState];
        double hVal = bestH[curState];
        double fVal = gVal + hVal;

        PathNode* node = createPathNode(adapter.formatName(curState, prevState), fVal, gVal, hVal, nullptr);
        if (head == nullptr) {
            head = node;
            tail = node;
        } else {
            tail->next = node;
            tail = node;
        }
    }

    return head;
}

// ============================================================================
// TASK 1 HELPERS: Mesh Network Packet Routing with Congestion Penalty
// ============================================================================

// A3: Degree of node v = number of outgoing edges (weightMatrix[v][j] > 0, j != v).
// Isolated here so it can be swapped for in-degree or symmetric degree in one place.
int degreeOf(double weightMatrix[100][100], int v) {
    int deg = 0;
    for (int j = 0; j < MAX_NODES; ++j) {
        if (j != v && hasEdge(weightMatrix[v][j])) {
            ++deg;
        }
    }
    return deg;
}

// R9: Isolated naming helper for Task 1.
std::string formatNameMesh(int node) {
    return std::to_string(node);
}

struct MeshAdapter {
    double (*weightMatrix)[100];
    int goalNode;
    int degrees[MAX_NODES];
    double minStep;

    bool isGoal(int state) const {
        return state == goalNode;
    }

    // Admissible & consistent heuristic: 0 at goal, minimum single-step cost elsewhere.
    double heuristic(int state) const {
        return (state == goalNode) ? 0.0 : minStep;
    }

    template <typename Callback>
    void forEachNeighbor(int u, Callback&& visit) const {
        for (int v = 0; v < MAX_NODES; ++v) {
            if (v != u && hasEdge(weightMatrix[u][v])) {
                // Congestion penalty applies to the target node v.
                double cost = weightMatrix[u][v] + 0.1 * static_cast<double>(degrees[v]);
                visit(v, cost);
            }
        }
    }

    std::string formatName(int state, int /*prevState*/) const {
        return formatNameMesh(state);
    }
};

// ============================================================================
// TASK 2 HELPERS: AUV Submarine Navigation in 3D Space
// ============================================================================

// R9: Isolated naming helper for Task 2.
std::string formatNameAUV(int node) {
    return std::to_string(node);
}

// 3D spatial heuristic for Task 2.
// Mode 1: Manhattan, Mode 2: Euclidean, Mode 3: Chebyshev. Invalid falls back to Mode 1.
double heuristicAUV(int coords[100][3], int node, int goal, int mode) {
    double dx = std::fabs(static_cast<double>(coords[node][0] - coords[goal][0]));
    double dy = std::fabs(static_cast<double>(coords[node][1] - coords[goal][1]));
    double dz = std::fabs(static_cast<double>(coords[node][2] - coords[goal][2]));

    if (mode == 2) {
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
    if (mode == 3) {
        return std::max(dx, std::max(dy, dz));
    }
    // Mode 1 (and fallback for any invalid mode): 3D Manhattan.
    return dx + dy + dz;
}

struct AUVAdapter {
    double (*weightMatrix)[100];
    int (*coords)[3];
    int goalPoint;
    int mode;

    bool isGoal(int state) const {
        return state == goalPoint;
    }

    double heuristic(int state) const {
        return heuristicAUV(coords, state, goalPoint, mode);
    }

    template <typename Callback>
    void forEachNeighbor(int u, Callback&& visit) const {
        for (int v = 0; v < MAX_NODES; ++v) {
            if (v != u && hasEdge(weightMatrix[u][v])) {
                // Vertical water-pressure depth penalty using floating-point division.
                double depthDiff = std::fabs(static_cast<double>(coords[v][2] - coords[u][2]));
                double cost = weightMatrix[u][v] * (1.0 + depthDiff / 100.0);
                visit(v, cost);
            }
        }
    }

    std::string formatName(int state, int /*prevState*/) const {
        return formatNameAUV(state);
    }
};

// ============================================================================
// TASK 3 HELPERS: Probabilistic Pathfinding under Hazard Risk
// ============================================================================

// R9: Isolated naming helper for Task 3.
std::string formatNameProbabilistic(int node) {
    return std::to_string(node);
}

// A5: 2D spatial heuristic for Task 3.
// Mode 1: Manhattan, Mode 2: Euclidean. Invalid falls back to Mode 1.
double heuristicProbabilistic(int coords[100][2], int node, int goal, int mode) {
    double dx = std::fabs(static_cast<double>(coords[node][0] - coords[goal][0]));
    double dy = std::fabs(static_cast<double>(coords[node][1] - coords[goal][1]));

    if (mode == 2) {
        return std::sqrt(dx * dx + dy * dy);
    }
    // Mode 1 (and fallback): 2D Manhattan.
    return dx + dy;
}

struct ProbabilisticAdapter {
    double (*distMatrix)[100];
    double (*riskMatrix)[100];
    int (*coords)[2];
    int goalPoint;
    double lambda;
    int mode;

    bool isGoal(int state) const {
        return state == goalPoint;
    }

    double heuristic(int state) const {
        return heuristicProbabilistic(coords, state, goalPoint, mode);
    }

    template <typename Callback>
    void forEachNeighbor(int u, Callback&& visit) const {
        for (int v = 0; v < MAX_NODES; ++v) {
            if (v != u && hasEdge(distMatrix[u][v])) {
                // Expected traversal cost combines distance and lambda-weighted hazard risk.
                double cost = distMatrix[u][v] + lambda * riskMatrix[u][v];
                visit(v, cost);
            }
        }
    }

    std::string formatName(int state, int /*prevState*/) const {
        return formatNameProbabilistic(state);
    }
};

// ============================================================================
// SHARED GRID CONVENTIONS (Tasks 4 & 5)
// ============================================================================

// A4: x is row (grid[x][y]), y is column; Up=x-1, Down=x+1, Left=y-1, Right=y+1.
// Indices 0..3 are cardinal (used by both Task 4 and Task 5); 4..7 are diagonal (Task 4).
const int DIR_DX[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
const int DIR_DY[8] = {0, 0, -1, 1, -1, 1, -1, 1};
const char* const DIR_NAMES[8] = {
    "Up", "Down", "Left", "Right",
    "Up-Left", "Up-Right", "Down-Left", "Down-Right"
};

inline bool inGridBounds(int x, int y, int m, int n) {
    return x >= 0 && x < m && y >= 0 && y < n;
}

inline int cellToId(int x, int y, int n) {
    return x * n + y;
}

// ============================================================================
// TASK 4 HELPERS: Game Character Navigation with Multi-Terrain Obstacles
// ============================================================================

// R9: Isolated naming helper for Task 4.
// First node is "Start"; every later node is named by the direction taken to reach it.
std::string formatNameGame(int curCell, int prevCell, int n) {
    if (prevCell < 0) {
        return "Start";
    }
    int px = prevCell / n;
    int py = prevCell % n;
    int cx = curCell / n;
    int cy = curCell % n;
    int dx = cx - px;
    int dy = cy - py;
    for (int d = 0; d < 8; ++d) {
        if (DIR_DX[d] == dx && DIR_DY[d] == dy) {
            return std::string(DIR_NAMES[d]);
        }
    }
    return "Start";
}

// Mode 1: Manhattan (not admissible with diagonal cost 1.5), Mode 2: Chebyshev (admissible).
double heuristicGrid2D(int x, int y, int gx, int gy, int mode) {
    double dx = std::fabs(static_cast<double>(x - gx));
    double dy = std::fabs(static_cast<double>(y - gy));
    if (mode == 2) {
        return std::max(dx, dy);
    }
    return dx + dy;
}

struct GameAdapter {
    int (*mapGrid)[100];
    int m;
    int n;
    int goalX;
    int goalY;
    int goalCell;
    int mode;

    bool isGoal(int state) const {
        return state == goalCell;
    }

    double heuristic(int state) const {
        int x = state / n;
        int y = state % n;
        return heuristicGrid2D(x, y, goalX, goalY, mode);
    }

    template <typename Callback>
    void forEachNeighbor(int state, Callback&& visit) const {
        int x = state / n;
        int y = state % n;

        for (int d = 0; d < 8; ++d) {
            int dx = DIR_DX[d];
            int dy = DIR_DY[d];
            int nx = x + dx;
            int ny = y + dy;

            if (!inGridBounds(nx, ny, m, n)) {
                continue;
            }
            // Value 1 is an impassable obstacle.
            if (mapGrid[nx][ny] == 1) {
                continue;
            }

            bool isDiagonal = (dx != 0 && dy != 0);
            if (isDiagonal) {
                // Corner-cutting rule: both adjacent cardinal cells must be in bounds and non-obstacles.
                if (!inGridBounds(x + dx, y, m, n) || !inGridBounds(x, y + dy, m, n)) {
                    continue;
                }
                if (mapGrid[x + dx][y] == 1 || mapGrid[x][y + dy] == 1) {
                    continue;
                }
            }

            // Terrain cost of the target cell: 0 -> flat (1.0), 2 -> swamp (2.0).
            double baseCost = (mapGrid[nx][ny] == 2) ? 2.0 : 1.0;
            double stepCost = isDiagonal ? (baseCost * 1.5) : baseCost;
            visit(cellToId(nx, ny, n), stepCost);
        }
    }

    std::string formatName(int state, int prevState) const {
        return formatNameGame(state, prevState, n);
    }
};

// ============================================================================
// TASK 5 HELPERS: VLSI PCB Wire Routing with Bending Penalty
// ============================================================================

const int PCB_NUM_DIRS = 5;    // 0..3 = Up/Down/Left/Right, 4 = None (start state)
const int PCB_DIR_NONE = 4;
const double PCB_BEND_PENALTY = 0.5;

// A6: In pcbGrid, 0 = free, 1 (non-zero) = obstacle.
inline bool isFreePCB(int cellValue) {
    return cellValue == 0;
}

// R9: Isolated naming helper for Task 5 -> "(x,y)".
std::string formatNamePCB(int cellId, int n) {
    int x = cellId / n;
    int y = cellId % n;
    return "(" + std::to_string(x) + "," + std::to_string(y) + ")";
}

struct PCBAdapter {
    double (*weightMatrix)[100];
    int m;
    int n;
    int exitX;
    int exitY;
    int exitCell;
    int mode;

    bool isGoal(int state) const {
        // Accept reaching the exit cell from any arrival direction.
        return (state / PCB_NUM_DIRS) == exitCell;
    }

    double heuristic(int state) const {
        int cell = state / PCB_NUM_DIRS;
        int x = cell / n;
        int y = cell % n;
        return heuristicGrid2D(x, y, exitX, exitY, mode);
    }

    template <typename Callback>
    void forEachNeighbor(int state, Callback&& visit) const {
        int cell = state / PCB_NUM_DIRS;
        int prevDir = state % PCB_NUM_DIRS;
        int x = cell / n;
        int y = cell % n;

        // 4-directional movement using the populated graph weightMatrix.
        for (int d = 0; d < 4; ++d) {
            int nx = x + DIR_DX[d];
            int ny = y + DIR_DY[d];
            if (!inGridBounds(nx, ny, m, n)) {
                continue;
            }
            int nextCell = cellToId(nx, ny, n);
            if (!hasEdge(weightMatrix[cell][nextCell])) {
                continue;
            }
            // Add bending penalty B = 0.5 when changing direction (first move has no penalty).
            double bend = (prevDir != PCB_DIR_NONE && prevDir != d) ? PCB_BEND_PENALTY : 0.0;
            double stepCost = weightMatrix[cell][nextCell] + bend;
            int nextState = nextCell * PCB_NUM_DIRS + d;
            visit(nextState, stepCost);
        }
    }

    std::string formatName(int state, int /*prevState*/) const {
        return formatNamePCB(state / PCB_NUM_DIRS, n);
    }
};

} // namespace

// R10: Single allowed global variable for Task 3 P_safe reporting.
double g_lastPsafe = 0.0;

// ============================================================================
// TASK 1: Mesh Network Packet Routing with Congestion Penalty
// ============================================================================
PathNode* findMeshPath(double weightMatrix[100][100], int startNode, int goalNode) {
    if (startNode < 0 || startNode >= MAX_NODES || goalNode < 0 || goalNode >= MAX_NODES) {
        return nullptr;
    }

    MeshAdapter adapter;
    adapter.weightMatrix = weightMatrix;
    adapter.goalNode = goalNode;

    // Precompute node degrees once (O(V^2)).
    for (int v = 0; v < MAX_NODES; ++v) {
        adapter.degrees[v] = degreeOf(weightMatrix, v);
    }

    // Precompute minimum possible single-edge cost across the entire graph.
    double minStep = INF;
    bool foundAnyEdge = false;
    for (int u = 0; u < MAX_NODES; ++u) {
        for (int v = 0; v < MAX_NODES; ++v) {
            if (u != v && hasEdge(weightMatrix[u][v])) {
                double stepCost = weightMatrix[u][v] + 0.1 * static_cast<double>(adapter.degrees[v]);
                if (stepCost < minStep) {
                    minStep = stepCost;
                }
                foundAnyEdge = true;
            }
        }
    }
    adapter.minStep = foundAnyEdge ? minStep : 0.0;

    return runAStar(MAX_NODES, startNode, adapter);
}

// ============================================================================
// TASK 2: AUV Submarine Navigation in 3D Space
// ============================================================================
PathNode* findAUVPath(double weightMatrix[100][100], int coords[100][3],
                      int startPoint, int goalPoint, int mode) {
    if (startPoint < 0 || startPoint >= MAX_NODES || goalPoint < 0 || goalPoint >= MAX_NODES) {
        return nullptr;
    }

    AUVAdapter adapter;
    adapter.weightMatrix = weightMatrix;
    adapter.coords = coords;
    adapter.goalPoint = goalPoint;
    adapter.mode = (mode == 1 || mode == 2 || mode == 3) ? mode : 1;

    return runAStar(MAX_NODES, startPoint, adapter);
}

// ============================================================================
// TASK 3: Probabilistic Pathfinding under Hazard Risk
// ============================================================================
PathNode* findProbabilisticPath(double distMatrix[100][100],
                                double riskMatrix[100][100], int coords[100][2],
                                int startPoint, int goalPoint, double lambda, int mode) {
    // Reset global survival probability at the start of every call.
    g_lastPsafe = 0.0;

    if (startPoint < 0 || startPoint >= MAX_NODES || goalPoint < 0 || goalPoint >= MAX_NODES) {
        return nullptr;
    }

    ProbabilisticAdapter adapter;
    adapter.distMatrix = distMatrix;
    adapter.riskMatrix = riskMatrix;
    adapter.coords = coords;
    adapter.goalPoint = goalPoint;
    adapter.lambda = lambda;
    adapter.mode = (mode == 1 || mode == 2) ? mode : 1;

    std::vector<int> pathStates;
    PathNode* path = runAStar(MAX_NODES, startPoint, adapter, &pathStates);
    if (path == nullptr) {
        g_lastPsafe = 0.0;
        return nullptr;
    }

    // Independent survival probability along the returned path (1.0 for single-node path).
    double pSafe = 1.0;
    for (std::size_t i = 1; i < pathStates.size(); ++i) {
        int u = pathStates[i - 1];
        int v = pathStates[i];
        pSafe *= (1.0 - riskMatrix[u][v]);
    }
    g_lastPsafe = pSafe;

    return path;
}

// ============================================================================
// TASK 4: Game Character Navigation with Multi-Terrain Obstacles
// ============================================================================
PathNode* findGamePath(int mapGrid[100][100], int m, int n, int startX,
                       int startY, int goalX, int goalY, int mode) {
    if (m <= 0 || m > 100 || n <= 0 || n > 100) {
        return nullptr;
    }
    if (!inGridBounds(startX, startY, m, n) || !inGridBounds(goalX, goalY, m, n)) {
        return nullptr;
    }
    // Start or goal on an obstacle is invalid (R6).
    if (mapGrid[startX][startY] == 1 || mapGrid[goalX][goalY] == 1) {
        return nullptr;
    }

    GameAdapter adapter;
    adapter.mapGrid = mapGrid;
    adapter.m = m;
    adapter.n = n;
    adapter.goalX = goalX;
    adapter.goalY = goalY;
    adapter.goalCell = cellToId(goalX, goalY, n);
    adapter.mode = (mode == 1 || mode == 2) ? mode : 1;

    return runAStar(m * n, cellToId(startX, startY, n), adapter);
}

// ============================================================================
// TASK 5: VLSI PCB Wire Routing with Bending Penalty
// ============================================================================
PathNode* findPCBPath(int pcbGrid[100][100], int m, int n, int startX,
                      int startY, int exitX, int exitY,
                      double weightMatrix[100][100], int mode) {
    // A6: weightMatrix is 100x100, so m * n must be in [1, 100].
    if (m <= 0 || n <= 0 || m * n > MAX_NODES) {
        return nullptr;
    }

    int numCells = m * n;

    // Step 1: Grid-to-graph transformation — zero the used m*n x m*n block first.
    for (int i = 0; i < numCells; ++i) {
        for (int j = 0; j < numCells; ++j) {
            weightMatrix[i][j] = 0.0;
        }
    }

    // Populate symmetric weight 1.0 between 4-adjacent free cells.
    for (int x = 0; x < m; ++x) {
        for (int y = 0; y < n; ++y) {
            if (!isFreePCB(pcbGrid[x][y])) {
                continue;
            }
            int u = cellToId(x, y, n);
            for (int d = 0; d < 4; ++d) {
                int nx = x + DIR_DX[d];
                int ny = y + DIR_DY[d];
                if (inGridBounds(nx, ny, m, n) && isFreePCB(pcbGrid[nx][ny])) {
                    int v = cellToId(nx, ny, n);
                    weightMatrix[u][v] = 1.0;
                }
            }
        }
    }

    // Step 2: Validate start and exit coordinates.
    if (!inGridBounds(startX, startY, m, n) || !inGridBounds(exitX, exitY, m, n)) {
        return nullptr;
    }
    if (!isFreePCB(pcbGrid[startX][startY]) || !isFreePCB(pcbGrid[exitX][exitY])) {
        return nullptr;
    }

    // Step 3: Run A* over the (cell, previous-direction) state space.
    PCBAdapter adapter;
    adapter.weightMatrix = weightMatrix;
    adapter.m = m;
    adapter.n = n;
    adapter.exitX = exitX;
    adapter.exitY = exitY;
    adapter.exitCell = cellToId(exitX, exitY, n);
    adapter.mode = (mode == 1 || mode == 2) ? mode : 1;

    int startState = cellToId(startX, startY, n) * PCB_NUM_DIRS + PCB_DIR_NONE;
    return runAStar(numCells * PCB_NUM_DIRS, startState, adapter);
}
