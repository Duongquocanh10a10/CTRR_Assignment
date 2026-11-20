#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

struct PathNode {
    std::string name;
    double f;
    double g;
    double h;
    PathNode* next;
};

PathNode* createPathNode(const std::string& name, double f, double g, double h, PathNode* next = nullptr) {
    PathNode* node = new PathNode;
    node->name = name;
    node->f = f;
    node->g = g;
    node->h = h;
    node->next = next;
    return node;
}

void freePath(PathNode* head) {
    while (head != nullptr) {
        PathNode* nextNode = head->next;
        delete head;
        head = nextNode;
    }
}

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

} // namespace

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

