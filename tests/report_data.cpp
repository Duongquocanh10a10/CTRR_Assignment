#include "../submission/Algo.h"
#include "../submission/PathNode.h"
#include "reference/ref_dijkstra.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

const double INF = 1e18;
const double EPS = 1e-9;

struct TraceEntry {
    double f, g, h;
    int state;
};

bool higherPrio(const TraceEntry& a, const TraceEntry& b) {
    if (std::fabs(a.f - b.f) > EPS) return a.f < b.f;
    if (std::fabs(a.h - b.h) > EPS) return a.h < b.h;
    return a.state < b.state;
}

// Instrumented A* matching submission/Algo.cpp's exact tie-breaking and lazy-deletion semantics,
// used ONLY to record expansion counts and step-by-step traces for the report.
template <typename NeighborFn, typename HeuristicFn, typename GoalFn>
int traceAStar(int numStates, int startState, GoalFn isGoal, HeuristicFn hFn, NeighborFn neighFn,
               bool printTrace, const std::function<std::string(int)>& labelFn) {
    std::vector<double> bestG(numStates, INF);
    std::vector<TraceEntry> heap;
    auto push = [&](const TraceEntry& e) {
        heap.push_back(e);
        std::size_t i = heap.size() - 1;
        while (i > 0) {
            std::size_t p = (i - 1) / 2;
            if (higherPrio(heap[i], heap[p])) { std::swap(heap[i], heap[p]); i = p; }
            else break;
        }
    };
    auto pop = [&]() {
        TraceEntry top = heap[0];
        heap[0] = heap.back();
        heap.pop_back();
        std::size_t i = 0, n = heap.size();
        while (true) {
            std::size_t l = 2 * i + 1, r = 2 * i + 2, b = i;
            if (l < n && higherPrio(heap[l], heap[b])) b = l;
            if (r < n && higherPrio(heap[r], heap[b])) b = r;
            if (b != i) { std::swap(heap[i], heap[b]); i = b; }
            else break;
        }
        return top;
    };

    double h0 = hFn(startState);
    bestG[startState] = 0.0;
    push({h0, 0.0, h0, startState});

    int expanded = 0;
    while (!heap.empty()) {
        TraceEntry cur = pop();
        if (cur.g > bestG[cur.state] + EPS) continue;
        ++expanded;
        if (printTrace) {
            std::cout << "  Step " << expanded << ": POP " << labelFn(cur.state)
                      << " (g=" << cur.g << ", h=" << cur.h << ", f=" << cur.f << ")\n";
        }
        if (isGoal(cur.state)) break;
        neighFn(cur.state, [&](int nxt, double cost) {
            double ng = cur.g + cost;
            if (ng < bestG[nxt] - EPS) {
                double nh = hFn(nxt);
                bestG[nxt] = ng;
                push({ng + nh, ng, nh, nxt});
                if (printTrace) {
                    std::cout << "    -> relax " << labelFn(nxt)
                              << " (stepCost=" << cost << ", newG=" << ng
                              << ", h=" << nh << ", f=" << (ng + nh) << ")\n";
                }
            }
        });
    }
    return expanded;
}

void printList(PathNode* head) {
    while (head != nullptr) {
        std::cout << "  Node: " << head->name
                  << " | g=" << head->g << " | h=" << head->h << " | f=" << head->f << "\n";
        head = head->next;
    }
}

void reportTask1() {
    std::cout << "=== REPORT DATA: TASK 1 (Mesh Trace) ===\n";
    double w[100][100];
    std::memset(w, 0, sizeof(w));
    // 5-node network (0..5), start=0, goal=4 (same as main.cpp Demo 1)
    w[0][1] = 1.0; w[0][2] = 1.4;
    w[1][2] = 1.0; w[1][3] = 1.0; w[1][4] = 2.0; w[1][5] = 1.0;
    w[2][4] = 1.6; w[3][4] = 1.5;

    int deg[100] = {0};
    for (int u = 0; u < 6; ++u) {
        for (int v = 0; v < 6; ++v) {
            if (u != v && w[u][v] > 0.0) ++deg[u];
        }
        std::cout << "  deg(" << u << ") = " << deg[u] << "\n";
    }
    double minStep = INF;
    for (int u = 0; u < 6; ++u)
        for (int v = 0; v < 6; ++v)
            if (u != v && w[u][v] > 0.0)
                minStep = std::min(minStep, w[u][v] + 0.1 * deg[v]);
    std::cout << "  minStep = " << minStep << "\n";

    traceAStar(
        6, 0,
        [](int s) { return s == 4; },
        [&](int s) { return (s == 4) ? 0.0 : minStep; },
        [&](int u, const std::function<void(int, double)>& visit) {
            for (int v = 0; v < 6; ++v)
                if (v != u && w[u][v] > 0.0)
                    visit(v, w[u][v] + 0.1 * deg[v]);
        },
        true,
        [](int s) { return std::to_string(s); });

    std::cout << "  Program findMeshPath output:\n";
    PathNode* p = findMeshPath(w, 0, 4);
    printList(p);
    freePath(p);
    std::cout << "\n";
}

void reportTask2() {
    std::cout << "=== REPORT DATA: TASK 2 (AUV 3D Heuristic Comparison) ===\n";
    // Test Case A: 50-node admissible 3D survey network (w >= 3D Manhattan)
    // Test Case B: 100-node admissible 3D survey network (w >= 3D Manhattan)
    // Test Case C: 50-node non-admissible 3D network (small arbitrary w in [1, 15], coords in [0, 50])
    struct CaseSpec {
        const char* name;
        int n;
        bool admissible;
        unsigned seed;
    };
    CaseSpec cases[3] = {
        {"Case A (50 nodes, admissible)", 50, true, 101},
        {"Case B (100 nodes, admissible)", 100, true, 202},
        {"Case C (50 nodes, non-admissible)", 50, false, 303}
    };

    double w[100][100];
    int c[100][3];

    for (const auto& cs : cases) {
        std::mt19937 rng(cs.seed);
        std::uniform_int_distribution<int> cd(0, 50);
        std::uniform_real_distribution<double> pr(0.0, 1.0);
        std::uniform_real_distribution<double> ex(0.5, 8.0);
        std::uniform_real_distribution<double> sw(1.0, 12.0);

        std::memset(w, 0, sizeof(w));
        std::memset(c, 0, sizeof(c));
        for (int i = 0; i < cs.n; ++i) {
            c[i][0] = cd(rng);
            c[i][1] = cd(rng);
            c[i][2] = cd(rng);
        }
        // Guarantee connectivity via a spine plus random edges
        for (int i = 0; i < cs.n; ++i) {
            for (int j = 0; j < cs.n; ++j) {
                if (i == j) continue;
                bool edge = (j == i + 1) || (pr(rng) < 0.12);
                if (edge) {
                    if (cs.admissible) {
                        double manh = std::fabs(double(c[i][0] - c[j][0])) +
                                      std::fabs(double(c[i][1] - c[j][1])) +
                                      std::fabs(double(c[i][2] - c[j][2]));
                        w[i][j] = std::max(1.0, manh) + ex(rng);
                    } else {
                        w[i][j] = sw(rng);
                    }
                }
            }
        }

        int s = 0, g = cs.n - 1;
        RefResult ref = refAUVDijkstra(w, c, s, g);
        std::cout << cs.name << " | Reference Dijkstra Optimal Cost = "
                  << std::fixed << std::setprecision(4) << ref.cost << "\n";

        const char* mNames[4] = {"", "Mode 1 (Manhattan)", "Mode 2 (Euclidean)", "Mode 3 (Chebyshev)"};
        for (int mode = 1; mode <= 3; ++mode) {
            // Measure average runtime over 200 runs of actual findAUVPath
            const int ITERS = 200;
            double finalCost = -1.0;
            auto t0 = std::chrono::steady_clock::now();
            for (int it = 0; it < ITERS; ++it) {
                PathNode* p = findAUVPath(w, c, s, g, mode);
                if (it == 0 && p != nullptr) {
                    PathNode* cur = p;
                    while (cur->next != nullptr) cur = cur->next;
                    finalCost = cur->g;
                }
                freePath(p);
            }
            auto t1 = std::chrono::steady_clock::now();
            double avgUs = std::chrono::duration<double, std::micro>(t1 - t0).count() / ITERS;

            // Count exact nodes expanded with matching comparator
            int expanded = traceAStar(
                cs.n, s,
                [&](int st) { return st == g; },
                [&](int st) {
                    double dx = std::fabs(double(c[st][0] - c[g][0]));
                    double dy = std::fabs(double(c[st][1] - c[g][1]));
                    double dz = std::fabs(double(c[st][2] - c[g][2]));
                    if (mode == 2) return std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (mode == 3) return std::max(dx, std::max(dy, dz));
                    return dx + dy + dz;
                },
                [&](int u, const std::function<void(int, double)>& visit) {
                    for (int v = 0; v < cs.n; ++v) {
                        if (v != u && w[u][v] > 0.0) {
                            double dz = std::fabs(double(c[v][2] - c[u][2]));
                            visit(v, w[u][v] * (1.0 + dz / 100.0));
                        }
                    }
                },
                false,
                [](int st) { return std::to_string(st); });

            double gap = finalCost - ref.cost;
            std::cout << "  " << mNames[mode]
                      << " | expanded=" << expanded
                      << " | avgTime=" << std::setprecision(2) << avgUs << " us"
                      << " | cost=" << std::setprecision(4) << finalCost
                      << " | gap=" << gap << "\n";
        }
    }
    std::cout << "\n";
}

void reportTask3() {
    std::cout << "=== REPORT DATA: TASK 3 (Hazard Risk Lambda Sweep) ===\n";
    double distM[100][100], riskM[100][100];
    int coords[100][2];
    std::memset(distM, 0, sizeof(distM));
    std::memset(riskM, 0, sizeof(riskM));
    std::memset(coords, 0, sizeof(coords));

    coords[0][0] = 0; coords[0][1] = 0;
    coords[1][0] = 4; coords[1][1] = 0;
    coords[2][0] = 4; coords[2][1] = 5;
    coords[3][0] = 8; coords[3][1] = 0;

    distM[0][1] = 5.0; riskM[0][1] = 0.40;
    distM[1][3] = 5.0; riskM[1][3] = 0.40;
    distM[0][2] = 9.0; riskM[0][2] = 0.05;
    distM[2][3] = 9.0; riskM[2][3] = 0.05;

    double lambdas[5] = {0.0, 5.0, 10.0, 12.0, 20.0};
    for (double lam : lambdas) {
        PathNode* p = findProbabilisticPath(distM, riskM, coords, 0, 3, lam, 2);
        std::string route;
        double totalG = 0.0;
        for (PathNode* cur = p; cur != nullptr; cur = cur->next) {
            if (!route.empty()) route += " -> ";
            route += cur->name;
            totalG = cur->g;
        }
        std::cout << "  lambda=" << std::setprecision(1) << lam
                  << " | route=" << route
                  << " | g=" << std::setprecision(4) << totalG
                  << " | P_safe=" << g_lastPsafe << "\n";
        freePath(p);
    }
    std::cout << "\n";
}

void reportTask4() {
    std::cout << "=== REPORT DATA: TASK 4 (Game Grid Trace, 4x4, Mode 2 Chebyshev) ===\n";
    int grid[100][100];
    std::memset(grid, 0, sizeof(grid));
    // 4x4 grid (same as main.cpp Demo 4)
    grid[0][1] = 1;
    grid[1][1] = 2;
    grid[1][2] = 2;
    grid[2][2] = 1;

    const int DX[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
    const int DY[8] = {0, 0, -1, 1, -1, 1, -1, 1};

    traceAStar(
        16, 0,
        [](int s) { return s == 15; },
        [](int s) {
            int x = s / 4, y = s % 4;
            return std::max(std::fabs(double(x - 3)), std::fabs(double(y - 3)));
        },
        [&](int s, const std::function<void(int, double)>& visit) {
            int x = s / 4, y = s % 4;
            for (int d = 0; d < 8; ++d) {
                int nx = x + DX[d], ny = y + DY[d];
                if (nx < 0 || nx >= 4 || ny < 0 || ny >= 4) continue;
                if (grid[nx][ny] == 1) continue;
                bool diag = (DX[d] != 0 && DY[d] != 0);
                if (diag && (grid[x + DX[d]][y] == 1 || grid[x][y + DY[d]] == 1)) {
                    std::cout << "    [corner-cut blocked] (" << x << "," << y << ") -> ("
                              << nx << "," << ny << ")\n";
                    continue;
                }
                double base = (grid[nx][ny] == 2) ? 2.0 : 1.0;
                visit(nx * 4 + ny, diag ? base * 1.5 : base);
            }
        },
        true,
        [](int s) { return "(" + std::to_string(s / 4) + "," + std::to_string(s % 4) + ")"; });

    std::cout << "  Program findGamePath (Mode 2) output:\n";
    PathNode* p = findGamePath(grid, 4, 4, 0, 0, 3, 3, 2);
    printList(p);
    freePath(p);
    std::cout << "\n";
}

} // namespace

int main() {
    reportTask1();
    reportTask2();
    reportTask3();
    reportTask4();
    return 0;
}
