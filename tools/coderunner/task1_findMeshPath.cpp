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

} // namespace

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

