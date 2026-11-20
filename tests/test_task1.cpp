#include "../submission/Algo.h"
#include "../submission/PathNode.h"
#include "reference/ref_dijkstra.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

const double EPS = 1e-7;
int g_testsPassed = 0;
int g_testsFailed = 0;

void check(bool cond, const std::string& msg) {
    if (cond) {
        ++g_testsPassed;
    } else {
        ++g_testsFailed;
        std::cerr << "[FAIL] " << msg << "\n";
    }
}

void zeroMatrix(double m[100][100]) {
    std::memset(m, 0, sizeof(double) * 100 * 100);
}

// Verify linked-list invariants: f == g + h, valid edge transitions, cost consistency.
bool verifyMeshPath(PathNode* head,
                    double weightMatrix[100][100],
                    int startNode,
                    int goalNode,
                    double& outTotalG) {
    if (head == nullptr) {
        return false;
    }
    int deg[100] = {0};
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            if (i != j && weightMatrix[i][j] > 0.0) {
                ++deg[i];
            }
        }
    }

    std::vector<PathNode*> nodes;
    for (PathNode* cur = head; cur != nullptr; cur = cur->next) {
        nodes.push_back(cur);
    }

    if (std::stoi(nodes.front()->name) != startNode) {
        return false;
    }
    if (std::stoi(nodes.back()->name) != goalNode) {
        return false;
    }
    if (std::fabs(nodes.front()->g) > EPS) {
        return false;
    }
    if (std::fabs(nodes.back()->h) > EPS) {
        return false;
    }

    double accumG = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (std::fabs(nodes[i]->f - (nodes[i]->g + nodes[i]->h)) > EPS) {
            return false;
        }
        if (i > 0) {
            int u = std::stoi(nodes[i - 1]->name);
            int v = std::stoi(nodes[i]->name);
            if (u == v || weightMatrix[u][v] <= 0.0) {
                return false;
            }
            accumG += weightMatrix[u][v] + 0.1 * static_cast<double>(deg[v]);
            if (std::fabs(nodes[i]->g - accumG) > EPS) {
                return false;
            }
        }
    }

    outTotalG = nodes.back()->g;
    return true;
}

void testEdgeCases() {
    double w[100][100];
    zeroMatrix(w);

    // 1. start == goal on empty graph
    {
        PathNode* p = findMeshPath(w, 5, 5);
        check(p != nullptr, "start == goal returns non-null");
        if (p != nullptr) {
            check(p->name == "5", "start == goal name == \"5\"");
            check(std::fabs(p->g) < EPS && std::fabs(p->h) < EPS && std::fabs(p->f) < EPS,
                  "start == goal has f=g=h=0");
            check(p->next == nullptr, "start == goal has single node");
            freePath(p);
        }
    }

    // 2. No path (disconnected)
    {
        w[0][1] = 2.0;
        w[2][3] = 3.0;
        PathNode* p = findMeshPath(w, 0, 3);
        check(p == nullptr, "disconnected nodes return nullptr");
    }

    // 3. Out of bounds
    {
        check(findMeshPath(w, -1, 3) == nullptr, "negative start returns nullptr");
        check(findMeshPath(w, 0, 100) == nullptr, "goal == 100 returns nullptr");
    }

    // 4. Self-loop ignored & zero weight means no edge
    {
        zeroMatrix(w);
        w[0][0] = 1.0; // self-loop should be ignored
        w[0][1] = 0.0; // zero weight = no edge
        PathNode* p = findMeshPath(w, 0, 1);
        check(p == nullptr, "zero weight is no edge");
    }

    // 5. Congestion penalty flips optimal route:
    // Direct path 0 -> 1 -> 3 has low base latency, but node 1 has high degree (e.g. 20),
    // adding +2.0 penalty so 0 -> 2 -> 3 wins!
    {
        zeroMatrix(w);
        w[0][1] = 1.0;
        w[1][3] = 1.0;
        // Give node 1 high out-degree (20 outgoing edges)
        for (int k = 10; k < 29; ++k) {
            w[1][k] = 5.0;
        }
        // Alternative path 0 -> 2 -> 3 with slightly higher base latency, low degree
        w[0][2] = 1.5;
        w[2][3] = 1.5;

        // Path 0->1->3 cost: (1.0 + 0.1*20) + (1.0 + 0.1*0) = 4.0
        // Path 0->2->3 cost: (1.5 + 0.1*1) + (1.5 + 0.1*0) = 3.1
        PathNode* p = findMeshPath(w, 0, 3);
        double totalG = -1.0;
        check(verifyMeshPath(p, w, 0, 3, totalG), "congestion flip path valid");
        check(std::fabs(totalG - 3.1) < EPS, "congestion flip chooses 0->2->3 with cost 3.1");
        if (p != nullptr && p->next != nullptr) {
            check(p->next->name == "2", "congestion flip middle node is 2");
        }
        freePath(p);
    }

    // 6. Tie-breaking test: 0->1->3 and 0->2->3 have identical costs and identical h.
    // Tie-breaker R8 chooses smaller node index (1 before 2).
    {
        zeroMatrix(w);
        w[0][1] = 2.0;
        w[1][3] = 2.0;
        w[0][2] = 2.0;
        w[2][3] = 2.0;

        PathNode* p = findMeshPath(w, 0, 3);
        double totalG = -1.0;
        check(verifyMeshPath(p, w, 0, 3, totalG), "tie-break path valid");
        if (p != nullptr && p->next != nullptr) {
            check(p->next->name == "1", "tie-break chooses smaller node index 1 over 2");
        }
        freePath(p);
    }

    // 7. Maximum 100-node chain (0 -> 1 -> ... -> 99)
    {
        zeroMatrix(w);
        for (int i = 0; i < 99; ++i) {
            w[i][i + 1] = 1.0;
        }
        PathNode* p = findMeshPath(w, 0, 99);
        RefResult ref = refMeshDijkstra(w, 0, 99);
        double totalG = -1.0;
        check(verifyMeshPath(p, w, 0, 99, totalG), "100-node chain path valid");
        check(ref.found && std::fabs(totalG - ref.cost) < EPS, "100-node chain matches reference Dijkstra");
        freePath(p);
    }
}

void testRandomAgainstReference() {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> probDist(0.0, 1.0);
    std::uniform_real_distribution<double> weightDist(0.5, 20.0);

    double w[100][100];
    int matched = 0;
    int totalInstances = 100;

    for (int tc = 0; tc < totalInstances; ++tc) {
        zeroMatrix(w);
        int n = (tc < 80) ? (5 + (tc % 20)) : 100;
        double edgeProb = (n <= 25) ? 0.25 : 0.08;

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j && probDist(rng) < edgeProb) {
                    w[i][j] = weightDist(rng);
                }
            }
        }

        int s = static_cast<int>(rng() % n);
        int g = static_cast<int>(rng() % n);

        PathNode* p = findMeshPath(w, s, g);
        RefResult ref = refMeshDijkstra(w, s, g);

        if (!ref.found) {
            check(p == nullptr, "random tc=" + std::to_string(tc) + " unreachable matches ref");
            if (p == nullptr) {
                ++matched;
            }
        } else {
            double totalG = -1.0;
            bool valid = verifyMeshPath(p, w, s, g, totalG);
            bool eq = valid && (std::fabs(totalG - ref.cost) < EPS);
            check(eq, "random tc=" + std::to_string(tc) + " optimal cost matches ref");
            if (eq) {
                ++matched;
            }
        }
        freePath(p);
    }

    std::cout << "Random Task 1 instances matched reference Dijkstra: "
              << matched << "/" << totalInstances << "\n";
}

} // namespace

int main() {
    testEdgeCases();
    testRandomAgainstReference();

    std::cout << "Task 1 Tests Summary: " << g_testsPassed << " passed, "
              << g_testsFailed << " failed.\n";
    return (g_testsFailed == 0) ? 0 : 1;
}
