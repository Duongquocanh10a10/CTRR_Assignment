#include "Algo.h"
#include "PathNode.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Demo copy of printPath (mirrors the LMS CodeRunner version; main.cpp is not sent to CodeRunner).
void printPath(PathNode* head) {
    std::cout << "Solution Path:\n";
    while (head != nullptr) {
        std::cout << "Node: " << head->name
                  << " | f: " << head->f
                  << " | g: " << head->g
                  << " | h: " << head->h << "\n";
        head = head->next;
    }
}

namespace {

void zeroDoubleMatrix(double m[100][100]) {
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            m[i][j] = 0.0;
        }
    }
}

void zeroIntGrid(int g[100][100]) {
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            g[i][j] = 0;
        }
    }
}

// ============================================================================
// Fig. 1 Example from Assignment Specification (S -> A -> D -> G, optimal cost = 4)
// ============================================================================
void demoFig1Example() {
    std::cout << "============================================================\n";
    std::cout << "DEMO 0: Fig. 1 A* Example (S -> A -> D -> G, optimal cost = 4)\n";
    std::cout << "============================================================\n";

    // Nodes: 0=S, 1=A, 2=B, 3=C, 4=D, 5=G
    const char* names[6] = {"S", "A", "B", "C", "D", "G"};
    const double hVals[6] = {3.0, 2.0, 2.0, 1.0, 1.0, 0.0};
    double adj[6][6];
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            adj[i][j] = 0.0;
        }
    }
    adj[0][1] = 1.0; // S -> A (1)
    adj[0][2] = 4.0; // S -> B (4)
    adj[1][3] = 2.0; // A -> C (2)
    adj[1][4] = 2.0; // A -> D (2)
    adj[2][4] = 1.0; // B -> D (1)
    adj[3][5] = 3.0; // C -> G (3)
    adj[4][5] = 1.0; // D -> G (1)

    const double INF = 1e18;
    const double EPS = 1e-9;
    std::vector<double> bestG(6, INF);
    std::vector<int> parent(6, -1);
    std::vector<bool> closed(6, false);

    bestG[0] = 0.0;
    while (true) {
        int bestU = -1;
        double bestF = INF;
        for (int i = 0; i < 6; ++i) {
            if (!closed[i] && bestG[i] < INF * 0.5) {
                double f = bestG[i] + hVals[i];
                if (bestU == -1 || f < bestF - EPS ||
                    (std::fabs(f - bestF) <= EPS && hVals[i] < hVals[bestU] - EPS)) {
                    bestU = i;
                    bestF = f;
                }
            }
        }
        if (bestU == -1 || bestU == 5) {
            break;
        }
        closed[bestU] = true;
        for (int v = 0; v < 6; ++v) {
            if (adj[bestU][v] > 0.0) {
                double ng = bestG[bestU] + adj[bestU][v];
                if (ng < bestG[v] - EPS) {
                    bestG[v] = ng;
                    parent[v] = bestU;
                    closed[v] = false;
                }
            }
        }
    }

    std::vector<int> order;
    for (int u = 5; u != -1; u = parent[u]) {
        order.push_back(u);
    }
    std::reverse(order.begin(), order.end());

    PathNode* head = nullptr;
    PathNode* tail = nullptr;
    for (std::size_t i = 0; i < order.size(); ++i) {
        int u = order[i];
        PathNode* node = createPathNode(names[u], bestG[u] + hVals[u], bestG[u], hVals[u], nullptr);
        if (head == nullptr) {
            head = tail = node;
        } else {
            tail->next = node;
            tail = node;
        }
    }
    printPath(head);
    freePath(head);
    std::cout << "\n";
}

// ============================================================================
// Demo 1: Mesh Network Packet Routing with Congestion Penalty
// ============================================================================
void demoTask1() {
    std::cout << "============================================================\n";
    std::cout << "DEMO 1: Task 1 - Mesh Network Routing (findMeshPath)\n";
    std::cout << "============================================================\n";

    double w[100][100];
    zeroDoubleMatrix(w);

    // 5-node mesh network (0 -> 4):
    // Node 1 is a congested hub (out-degree 4: edges to 2, 3, 4, 5).
    // Path 0 -> 1 -> 4 has base latency 1.0 + 2.0 = 3.0, plus congestion 0.1*4 = 0.4 => 3.4.
    // Path 0 -> 2 -> 4 has base latency 1.4 + 1.6 = 3.0, plus congestion 0.1*1 = 0.1 => 3.1.
    w[0][1] = 1.0;
    w[0][2] = 1.4;
    w[1][2] = 1.0;
    w[1][3] = 1.0;
    w[1][4] = 2.0;
    w[1][5] = 1.0;
    w[2][4] = 1.6;
    w[3][4] = 1.5;

    PathNode* path = findMeshPath(w, 0, 4);
    printPath(path);
    freePath(path);
    std::cout << "\n";
}

// ============================================================================
// Demo 2: AUV Submarine Navigation in 3D Space
// ============================================================================
void demoTask2() {
    std::cout << "============================================================\n";
    std::cout << "DEMO 2: Task 2 - AUV 3D Navigation (findAUVPath)\n";
    std::cout << "============================================================\n";

    double w[100][100];
    zeroDoubleMatrix(w);
    int coords[100][3];
    for (int i = 0; i < 100; ++i) {
        coords[i][0] = coords[i][1] = coords[i][2] = 0;
    }

    // 4 stations: 0(0,0,0), 1(5,0,200) [deep], 2(5,2,10) [shallow], 3(10,0,0)
    coords[0][0] = 0;  coords[0][1] = 0; coords[0][2] = 0;
    coords[1][0] = 5;  coords[1][1] = 0; coords[1][2] = 200;
    coords[2][0] = 5;  coords[2][1] = 2; coords[2][2] = 10;
    coords[3][0] = 10; coords[3][1] = 0; coords[3][2] = 0;

    w[0][1] = 10.0; w[1][3] = 10.0; // Deep dive: 10*(1+2) + 10*(1+2) = 60.0
    w[0][2] = 12.0; w[2][3] = 12.0; // Shallow:   12*(1+0.1) + 12*(1+0.1) = 26.4

    const char* modeNames[4] = {"", "Mode 1 (3D Manhattan)", "Mode 2 (3D Euclidean)", "Mode 3 (3D Chebyshev)"};
    for (int mode = 1; mode <= 3; ++mode) {
        std::cout << "-- " << modeNames[mode] << " --\n";
        PathNode* path = findAUVPath(w, coords, 0, 3, mode);
        printPath(path);
        freePath(path);
    }
    std::cout << "\n";
}

// ============================================================================
// Demo 3: Probabilistic Pathfinding under Hazard Risk
// ============================================================================
void demoTask3() {
    std::cout << "============================================================\n";
    std::cout << "DEMO 3: Task 3 - Probabilistic Pathfinding (findProbabilisticPath)\n";
    std::cout << "============================================================\n";

    double distM[100][100];
    double riskM[100][100];
    zeroDoubleMatrix(distM);
    zeroDoubleMatrix(riskM);

    int coords[100][2];
    for (int i = 0; i < 100; ++i) {
        coords[i][0] = coords[i][1] = 0;
    }

    // Route A (short, risky): 0 -> 1 -> 3 (dist=10, risk=0.40 each, P_safe=0.36)
    // Route B (long, safe):   0 -> 2 -> 3 (dist=18, risk=0.05 each, P_safe=0.9025)
    coords[0][0] = 0; coords[0][1] = 0;
    coords[1][0] = 4; coords[1][1] = 0;
    coords[2][0] = 4; coords[2][1] = 5;
    coords[3][0] = 8; coords[3][1] = 0;

    distM[0][1] = 5.0; riskM[0][1] = 0.40;
    distM[1][3] = 5.0; riskM[1][3] = 0.40;
    distM[0][2] = 9.0; riskM[0][2] = 0.05;
    distM[2][3] = 9.0; riskM[2][3] = 0.05;

    double lambdas[3] = {0.0, 10.0, 20.0};
    for (int i = 0; i < 3; ++i) {
        std::cout << "-- lambda = " << lambdas[i] << " (Mode 2: Euclidean) --\n";
        PathNode* path = findProbabilisticPath(distM, riskM, coords, 0, 3, lambdas[i], 2);
        printPath(path);
        std::cout << "P_safe = " << g_lastPsafe << "\n";
        freePath(path);
    }
    std::cout << "\n";
}

// ============================================================================
// Demo 4: Game Character Navigation with Multi-Terrain Obstacles
// ============================================================================
void demoTask4() {
    std::cout << "============================================================\n";
    std::cout << "DEMO 4: Task 4 - Multi-Terrain Game Grid (findGamePath)\n";
    std::cout << "============================================================\n";

    int grid[100][100];
    zeroIntGrid(grid);
    // 4x4 grid: 0=flat, 1=obstacle, 2=swamp
    // Row 0: 0 1 0 0
    // Row 1: 0 2 2 0
    // Row 2: 0 0 1 0
    // Row 3: 0 0 0 0
    grid[0][1] = 1;
    grid[1][1] = 2;
    grid[1][2] = 2;
    grid[2][2] = 1;

    std::cout << "-- Mode 2 (Chebyshev, admissible) from (0,0) to (3,3) --\n";
    PathNode* path = findGamePath(grid, 4, 4, 0, 0, 3, 3, 2);
    printPath(path);
    freePath(path);
    std::cout << "\n";
}

// ============================================================================
// Demo 5: VLSI PCB Wire Routing with Bending Penalty
// ============================================================================
void demoTask5() {
    std::cout << "============================================================\n";
    std::cout << "DEMO 5: Task 5 - VLSI PCB Wire Routing (findPCBPath)\n";
    std::cout << "============================================================\n";

    int grid[100][100];
    zeroIntGrid(grid);
    for (int r = 0; r < 6; ++r) {
        for (int c = 0; c < 5; ++c) {
            grid[r][c] = 1;
        }
    }
    // Corridor A (8 steps, 7 bends = 11.5) vs Corridor B (10 steps, 2 bends = 11.0)
    int corrA[9][2] = {{0,0},{0,1},{1,1},{1,2},{2,2},{2,3},{3,3},{3,4},{4,4}};
    for (int i = 0; i < 9; ++i) {
        grid[corrA[i][0]][corrA[i][1]] = 0;
    }
    for (int r = 0; r <= 5; ++r) grid[r][0] = 0;
    for (int c = 0; c <= 4; ++c) grid[5][c] = 0;
    grid[4][4] = 0;

    double w[100][100];
    zeroDoubleMatrix(w);

    std::cout << "-- Mode 1 (Manhattan) from (0,0) to (4,4) on 6x5 PCB --\n";
    PathNode* path = findPCBPath(grid, 6, 5, 0, 0, 4, 4, w, 1);
    printPath(path);

    int edgeCount = 0;
    for (int i = 0; i < 30; ++i) {
        for (int j = 0; j < 30; ++j) {
            if (w[i][j] > 0.0) {
                ++edgeCount;
            }
        }
    }
    std::cout << "Populated weightMatrix (30x30) directed edges: " << edgeCount << "\n";
    freePath(path);
    std::cout << "\n";
}

} // namespace

int main() {
    demoFig1Example();
    demoTask1();
    demoTask2();
    demoTask3();
    demoTask4();
    demoTask5();
    return 0;
}
