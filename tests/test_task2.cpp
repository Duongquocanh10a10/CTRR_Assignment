#include "../submission/Algo.h"
#include "../submission/PathNode.h"
#include "reference/ref_dijkstra.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
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

void zeroCoords(int c[100][3]) {
    std::memset(c, 0, sizeof(int) * 100 * 3);
}

double expectedH(int coords[100][3], int u, int goal, int mode) {
    double dx = std::fabs(static_cast<double>(coords[u][0] - coords[goal][0]));
    double dy = std::fabs(static_cast<double>(coords[u][1] - coords[goal][1]));
    double dz = std::fabs(static_cast<double>(coords[u][2] - coords[goal][2]));
    if (mode == 2) return std::sqrt(dx * dx + dy * dy + dz * dz);
    if (mode == 3) return std::max(dx, std::max(dy, dz));
    return dx + dy + dz;
}

bool verifyAUVPath(PathNode* head,
                   double w[100][100],
                   int coords[100][3],
                   int startPoint,
                   int goalPoint,
                   int mode,
                   double& outTotalG) {
    if (head == nullptr) return false;
    int effMode = (mode == 1 || mode == 2 || mode == 3) ? mode : 1;

    std::vector<PathNode*> nodes;
    for (PathNode* cur = head; cur != nullptr; cur = cur->next) {
        nodes.push_back(cur);
    }

    if (std::stoi(nodes.front()->name) != startPoint) return false;
    if (std::stoi(nodes.back()->name) != goalPoint) return false;
    if (std::fabs(nodes.front()->g) > EPS) return false;
    if (std::fabs(nodes.back()->h) > EPS) return false;

    double accumG = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        int v = std::stoi(nodes[i]->name);
        double hExp = expectedH(coords, v, goalPoint, effMode);
        if (std::fabs(nodes[i]->h - hExp) > EPS) return false;
        if (std::fabs(nodes[i]->f - (nodes[i]->g + nodes[i]->h)) > EPS) return false;

        if (i > 0) {
            int u = std::stoi(nodes[i - 1]->name);
            if (u == v || w[u][v] <= 0.0) return false;
            double dz = std::fabs(static_cast<double>(coords[v][2] - coords[u][2]));
            accumG += w[u][v] * (1.0 + dz / 100.0);
            if (std::fabs(nodes[i]->g - accumG) > EPS) return false;
        }
    }

    outTotalG = nodes.back()->g;
    return true;
}

void testEdgeCases() {
    double w[100][100];
    int c[100][3];
    zeroMatrix(w);
    zeroCoords(c);

    // 1. start == goal
    {
        c[4][0] = 10; c[4][1] = 20; c[4][2] = 50;
        PathNode* p = findAUVPath(w, c, 4, 4, 2);
        check(p != nullptr, "start == goal returns non-null");
        if (p != nullptr) {
            check(p->name == "4" && std::fabs(p->g) < EPS && std::fabs(p->h) < EPS &&
                      std::fabs(p->f) < EPS && p->next == nullptr,
                  "start == goal single node with f=g=h=0");
            freePath(p);
        }
    }

    // 2. Disconnected & out-of-bounds
    {
        w[0][1] = 5.0;
        check(findAUVPath(w, c, 0, 2, 1) == nullptr, "disconnected returns nullptr");
        check(findAUVPath(w, c, -1, 1, 1) == nullptr, "negative start returns nullptr");
        check(findAUVPath(w, c, 0, 100, 1) == nullptr, "goal == 100 returns nullptr");
    }

    // 3. Invalid mode falls back to Mode 1 (3D Manhattan)
    {
        zeroMatrix(w);
        zeroCoords(c);
        c[0][0] = 0; c[0][1] = 0; c[0][2] = 0;
        c[1][0] = 3; c[1][1] = 4; c[1][2] = 12;
        w[0][1] = 20.0;

        PathNode* p = findAUVPath(w, c, 0, 1, 99);
        double totalG = -1.0;
        check(verifyAUVPath(p, w, c, 0, 1, 99, totalG), "invalid mode falls back to Mode 1");
        if (p != nullptr) {
            check(std::fabs(p->h - 19.0) < EPS, "invalid mode h(0) == 3+4+12 == 19.0");
            freePath(p);
        }
    }

    // 4. Depth penalty flips route:
    // 0 -> 1 -> 3 dives deep (z=0 -> z=200 -> z=0), penalty multiplies cost by 3.0 each leg.
    // 0 -> 2 -> 3 stays flat at z=0.
    {
        zeroMatrix(w);
        zeroCoords(c);
        c[0][0] = 0;  c[0][1] = 0;  c[0][2] = 0;
        c[1][0] = 5;  c[1][1] = 0;  c[1][2] = 200;
        c[2][0] = 5;  c[2][1] = 0;  c[2][2] = 0;
        c[3][0] = 10; c[3][1] = 0;  c[3][2] = 0;

        w[0][1] = 10.0; w[1][3] = 10.0; // cost = 10*(1+2) + 10*(1+2) = 60.0
        w[0][2] = 15.0; w[2][3] = 15.0; // cost = 15*1 + 15*1 = 30.0

        for (int mode = 1; mode <= 3; ++mode) {
            PathNode* p = findAUVPath(w, c, 0, 3, mode);
            double totalG = -1.0;
            check(verifyAUVPath(p, w, c, 0, 3, mode, totalG),
                  "depth-penalty flip valid in mode " + std::to_string(mode));
            check(std::fabs(totalG - 30.0) < EPS,
                  "depth-penalty flip picks flat route (30.0) in mode " + std::to_string(mode));
            freePath(p);
        }
    }

    // 5. Max 100-node graph
    {
        zeroMatrix(w);
        zeroCoords(c);
        for (int i = 0; i < 100; ++i) {
            c[i][0] = i;
            c[i][1] = (i % 5);
            c[i][2] = (i * 3) % 40;
            if (i < 99) {
                w[i][i + 1] = 50.0; // >= any 3D metric between adjacent nodes
            }
        }
        RefResult ref = refAUVDijkstra(w, c, 0, 99);
        for (int mode = 1; mode <= 3; ++mode) {
            PathNode* p = findAUVPath(w, c, 0, 99, mode);
            double totalG = -1.0;
            check(verifyAUVPath(p, w, c, 0, 99, mode, totalG),
                  "100-node AUV valid in mode " + std::to_string(mode));
            check(ref.found && std::fabs(totalG - ref.cost) < EPS,
                  "100-node AUV optimal in mode " + std::to_string(mode));
            freePath(p);
        }
    }
}

void testRandomAdmissibleAndNonAdmissible() {
    std::mt19937 rng(123);
    std::uniform_int_distribution<int> coordDist(0, 50);
    std::uniform_real_distribution<double> probDist(0.0, 1.0);
    std::uniform_real_distribution<double> extraDist(0.0, 10.0);
    std::uniform_real_distribution<double> smallWDist(1.0, 15.0);

    double w[100][100];
    int c[100][3];

    // Part A: Admissible instances (w(u,v) >= 3D Manhattan(u,v), so all 3 modes are admissible)
    int admissibleTotal = 60;
    int admissibleMatched[4] = {0, 0, 0, 0};

    for (int tc = 0; tc < admissibleTotal; ++tc) {
        zeroMatrix(w);
        zeroCoords(c);
        int n = 15 + (tc % 20);
        for (int i = 0; i < n; ++i) {
            c[i][0] = coordDist(rng);
            c[i][1] = coordDist(rng);
            c[i][2] = coordDist(rng);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j && probDist(rng) < 0.25) {
                    double manh = expectedH(c, i, j, 1);
                    w[i][j] = std::max(1.0, manh) + extraDist(rng);
                }
            }
        }
        int s = static_cast<int>(rng() % n);
        int g = static_cast<int>(rng() % n);
        RefResult ref = refAUVDijkstra(w, c, s, g);

        for (int mode = 1; mode <= 3; ++mode) {
            PathNode* p = findAUVPath(w, c, s, g, mode);
            if (!ref.found) {
                check(p == nullptr, "admissible AUV unreachable mode " + std::to_string(mode));
                if (p == nullptr) ++admissibleMatched[mode];
            } else {
                double totalG = -1.0;
                bool ok = verifyAUVPath(p, w, c, s, g, mode, totalG) &&
                          (std::fabs(totalG - ref.cost) < EPS);
                check(ok, "admissible AUV optimal tc=" + std::to_string(tc) +
                              " mode=" + std::to_string(mode));
                if (ok) ++admissibleMatched[mode];
            }
            freePath(p);
        }
    }

    std::cout << "Admissible AUV instances matched reference Dijkstra: "
              << "Mode1=" << admissibleMatched[1] << "/" << admissibleTotal
              << ", Mode2=" << admissibleMatched[2] << "/" << admissibleTotal
              << ", Mode3=" << admissibleMatched[3] << "/" << admissibleTotal << "\n";

    // Part B: General (potentially non-admissible) instances where w(u,v) can be smaller than h.
    // Do NOT assert equality; verify path validity and record optimality gaps.
    int generalTotal = 60;
    int suboptCount[4] = {0, 0, 0, 0};
    double maxGap[4] = {0.0, 0.0, 0.0, 0.0};

    for (int tc = 0; tc < generalTotal; ++tc) {
        zeroMatrix(w);
        zeroCoords(c);
        int n = 15 + (tc % 20);
        for (int i = 0; i < n; ++i) {
            c[i][0] = coordDist(rng);
            c[i][1] = coordDist(rng);
            c[i][2] = coordDist(rng);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j && probDist(rng) < 0.25) {
                    w[i][j] = smallWDist(rng);
                }
            }
        }
        int s = static_cast<int>(rng() % n);
        int g = static_cast<int>(rng() % n);
        RefResult ref = refAUVDijkstra(w, c, s, g);

        for (int mode = 1; mode <= 3; ++mode) {
            PathNode* p = findAUVPath(w, c, s, g, mode);
            if (!ref.found) {
                check(p == nullptr, "general AUV unreachable");
            } else {
                double totalG = -1.0;
                bool valid = verifyAUVPath(p, w, c, s, g, mode, totalG);
                check(valid, "general AUV path valid");
                double gap = totalG - ref.cost;
                check(gap >= -EPS, "general AUV cannot beat reference Dijkstra");
                if (gap > EPS) {
                    ++suboptCount[mode];
                    if (gap > maxGap[mode]) maxGap[mode] = gap;
                }
            }
            freePath(p);
        }
    }

    std::cout << std::fixed << std::setprecision(4)
              << "General (non-admissible) AUV gaps over " << generalTotal << " instances:\n"
              << "  Mode 1 (Manhattan):  " << suboptCount[1] << " suboptimal, maxGap=" << maxGap[1] << "\n"
              << "  Mode 2 (Euclidean):  " << suboptCount[2] << " suboptimal, maxGap=" << maxGap[2] << "\n"
              << "  Mode 3 (Chebyshev):  " << suboptCount[3] << " suboptimal, maxGap=" << maxGap[3] << "\n";
}

} // namespace

int main() {
    testEdgeCases();
    testRandomAdmissibleAndNonAdmissible();

    std::cout << "Task 2 Tests Summary: " << g_testsPassed << " passed, "
              << g_testsFailed << " failed.\n";
    return (g_testsFailed == 0) ? 0 : 1;
}
