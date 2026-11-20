#include "../submission/Algo.h"
#include "../submission/PathNode.h"
#include "reference/ref_dijkstra.h"

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

void zeroCoords(int c[100][2]) {
    std::memset(c, 0, sizeof(int) * 100 * 2);
}

double expectedH(int coords[100][2], int u, int goal, int mode) {
    double dx = std::fabs(static_cast<double>(coords[u][0] - coords[goal][0]));
    double dy = std::fabs(static_cast<double>(coords[u][1] - coords[goal][1]));
    if (mode == 2) return std::sqrt(dx * dx + dy * dy);
    return dx + dy;
}

bool verifyProbPath(PathNode* head,
                    double distM[100][100],
                    double riskM[100][100],
                    int coords[100][2],
                    int startPoint,
                    int goalPoint,
                    double lambda,
                    int mode,
                    double& outTotalG,
                    double& outExpectedPsafe) {
    if (head == nullptr) return false;
    int effMode = (mode == 1 || mode == 2) ? mode : 1;

    std::vector<PathNode*> nodes;
    for (PathNode* cur = head; cur != nullptr; cur = cur->next) {
        nodes.push_back(cur);
    }

    if (std::stoi(nodes.front()->name) != startPoint) return false;
    if (std::stoi(nodes.back()->name) != goalPoint) return false;
    if (std::fabs(nodes.front()->g) > EPS) return false;
    if (std::fabs(nodes.back()->h) > EPS) return false;

    double accumG = 0.0;
    double pSafe = 1.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        int v = std::stoi(nodes[i]->name);
        double hExp = expectedH(coords, v, goalPoint, effMode);
        if (std::fabs(nodes[i]->h - hExp) > EPS) return false;
        if (std::fabs(nodes[i]->f - (nodes[i]->g + nodes[i]->h)) > EPS) return false;

        if (i > 0) {
            int u = std::stoi(nodes[i - 1]->name);
            if (u == v || distM[u][v] <= 0.0) return false;
            accumG += distM[u][v] + lambda * riskM[u][v];
            pSafe *= (1.0 - riskM[u][v]);
            if (std::fabs(nodes[i]->g - accumG) > EPS) return false;
        }
    }

    outTotalG = nodes.back()->g;
    outExpectedPsafe = pSafe;
    return true;
}

void testEdgeCasesAndLambdaTradeoff() {
    double distM[100][100];
    double riskM[100][100];
    int coords[100][2];
    zeroMatrix(distM);
    zeroMatrix(riskM);
    zeroCoords(coords);

    // 1. start == goal -> single node, g_lastPsafe == 1.0
    {
        PathNode* p = findProbabilisticPath(distM, riskM, coords, 7, 7, 10.0, 2);
        check(p != nullptr, "start == goal returns non-null");
        check(std::fabs(g_lastPsafe - 1.0) < EPS, "start == goal sets g_lastPsafe = 1.0");
        if (p != nullptr) {
            check(p->name == "7" && std::fabs(p->g) < EPS && std::fabs(p->h) < EPS &&
                      std::fabs(p->f) < EPS && p->next == nullptr,
                  "start == goal single node with f=g=h=0");
            freePath(p);
        }
    }

    // 2. No path -> nullptr and g_lastPsafe reset to 0.0
    {
        distM[0][1] = 5.0;
        riskM[0][1] = 0.2;
        PathNode* p = findProbabilisticPath(distM, riskM, coords, 0, 3, 5.0, 1);
        check(p == nullptr, "disconnected returns nullptr");
        check(std::fabs(g_lastPsafe) < EPS, "disconnected resets g_lastPsafe = 0.0");
    }

    // 3. Out of bounds -> nullptr and g_lastPsafe = 0.0
    {
        g_lastPsafe = 0.99;
        check(findProbabilisticPath(distM, riskM, coords, -1, 1, 5.0, 1) == nullptr,
              "negative start returns nullptr");
        check(std::fabs(g_lastPsafe) < EPS, "out-of-bounds resets g_lastPsafe = 0.0");
    }

    // 4. Risk vs Length trade-off across lambda values (including lambda = 0)
    // Route A (short, risky): 0 -> 1 -> 3, dist = 5 + 5 = 10, risk = 0.4, 0.4
    //   g_A(lambda) = 10 + 0.8 * lambda, P_safe(A) = 0.6 * 0.6 = 0.36
    // Route B (long, safe):   0 -> 2 -> 3, dist = 9 + 9 = 18, risk = 0.05, 0.05
    //   g_B(lambda) = 18 + 0.1 * lambda, P_safe(B) = 0.95 * 0.95 = 0.9025
    // Threshold: 10 + 0.8*lambda = 18 + 0.1*lambda => 0.7*lambda = 8 => lambda ~= 11.4286
    {
        zeroMatrix(distM);
        zeroMatrix(riskM);
        zeroCoords(coords);
        coords[0][0] = 0;  coords[0][1] = 0;
        coords[1][0] = 4;  coords[1][1] = 0;
        coords[2][0] = 4;  coords[2][1] = 5;
        coords[3][0] = 8;  coords[3][1] = 0;

        distM[0][1] = 5.0;  riskM[0][1] = 0.40;
        distM[1][3] = 5.0;  riskM[1][3] = 0.40;
        distM[0][2] = 9.0;  riskM[0][2] = 0.05;
        distM[2][3] = 9.0;  riskM[2][3] = 0.05;

        // lambda = 0 -> Route A (cost 10.0, P_safe = 0.36)
        {
            PathNode* p = findProbabilisticPath(distM, riskM, coords, 0, 3, 0.0, 2);
            double totalG = -1.0, expPsafe = -1.0;
            check(verifyProbPath(p, distM, riskM, coords, 0, 3, 0.0, 2, totalG, expPsafe),
                  "lambda=0 path valid");
            check(std::fabs(totalG - 10.0) < EPS, "lambda=0 chooses short risky route (10.0)");
            check(std::fabs(g_lastPsafe - 0.36) < EPS, "lambda=0 g_lastPsafe == 0.36");
            freePath(p);
        }

        // lambda = 10 -> Route A still cheaper (18.0 vs 19.0)
        {
            PathNode* p = findProbabilisticPath(distM, riskM, coords, 0, 3, 10.0, 2);
            double totalG = -1.0, expPsafe = -1.0;
            check(verifyProbPath(p, distM, riskM, coords, 0, 3, 10.0, 2, totalG, expPsafe),
                  "lambda=10 path valid");
            check(std::fabs(totalG - 18.0) < EPS, "lambda=10 chooses Route A (18.0)");
            check(std::fabs(g_lastPsafe - 0.36) < EPS, "lambda=10 g_lastPsafe == 0.36");
            freePath(p);
        }

        // lambda = 20 -> Route B cheaper (20.0 vs 26.0, P_safe = 0.9025)
        {
            PathNode* p = findProbabilisticPath(distM, riskM, coords, 0, 3, 20.0, 2);
            double totalG = -1.0, expPsafe = -1.0;
            check(verifyProbPath(p, distM, riskM, coords, 0, 3, 20.0, 2, totalG, expPsafe),
                  "lambda=20 path valid");
            check(std::fabs(totalG - 20.0) < EPS, "lambda=20 switches to long safe route (20.0)");
            check(std::fabs(g_lastPsafe - 0.9025) < EPS, "lambda=20 g_lastPsafe == 0.9025");
            freePath(p);
        }
    }

    // 5. 100-node chain
    {
        zeroMatrix(distM);
        zeroMatrix(riskM);
        zeroCoords(coords);
        for (int i = 0; i < 100; ++i) {
            coords[i][0] = i;
            coords[i][1] = 0;
            if (i < 99) {
                distM[i][i + 1] = 2.0;
                riskM[i][i + 1] = 0.01;
            }
        }
        RefResult ref = refProbabilisticDijkstra(distM, riskM, 0, 99, 5.0);
        PathNode* p = findProbabilisticPath(distM, riskM, coords, 0, 99, 5.0, 2);
        double totalG = -1.0, expPsafe = -1.0;
        check(verifyProbPath(p, distM, riskM, coords, 0, 99, 5.0, 2, totalG, expPsafe),
              "100-node probabilistic path valid");
        check(ref.found && std::fabs(totalG - ref.cost) < EPS, "100-node matches reference");
        check(std::fabs(g_lastPsafe - std::pow(0.99, 99)) < EPS, "100-node g_lastPsafe == 0.99^99");
        freePath(p);
    }
}

void testRandomAgainstReference() {
    std::mt19937 rng(777);
    std::uniform_int_distribution<int> coordDist(0, 50);
    std::uniform_real_distribution<double> probDist(0.0, 1.0);
    std::uniform_real_distribution<double> riskDist(0.0, 0.95);
    std::uniform_real_distribution<double> extraDist(0.0, 10.0);
    std::uniform_real_distribution<double> smallDDist(1.0, 15.0);

    double distM[100][100];
    double riskM[100][100];
    int coords[100][2];

    // Part A: Admissible instances (distMatrix[u][v] >= 2D Manhattan(u,v), lambda >= 0)
    int admissibleTotal = 60;
    int matched[3] = {0, 0, 0};

    for (int tc = 0; tc < admissibleTotal; ++tc) {
        zeroMatrix(distM);
        zeroMatrix(riskM);
        zeroCoords(coords);
        int n = 15 + (tc % 20);
        double lambda = (tc % 4) * 5.0; // tests 0, 5, 10, 15

        for (int i = 0; i < n; ++i) {
            coords[i][0] = coordDist(rng);
            coords[i][1] = coordDist(rng);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j && probDist(rng) < 0.25) {
                    double manh = expectedH(coords, i, j, 1);
                    distM[i][j] = std::max(1.0, manh) + extraDist(rng);
                    riskM[i][j] = riskDist(rng);
                }
            }
        }

        int s = static_cast<int>(rng() % n);
        int g = static_cast<int>(rng() % n);
        RefResult ref = refProbabilisticDijkstra(distM, riskM, s, g, lambda);

        for (int mode = 1; mode <= 2; ++mode) {
            PathNode* p = findProbabilisticPath(distM, riskM, coords, s, g, lambda, mode);
            if (!ref.found) {
                check(p == nullptr && std::fabs(g_lastPsafe) < EPS,
                      "admissible Prob unreachable mode " + std::to_string(mode));
                if (p == nullptr) ++matched[mode];
            } else {
                double totalG = -1.0, expPsafe = -1.0;
                bool ok = verifyProbPath(p, distM, riskM, coords, s, g, lambda, mode, totalG, expPsafe) &&
                          (std::fabs(totalG - ref.cost) < EPS) &&
                          (std::fabs(g_lastPsafe - expPsafe) < EPS);
                check(ok, "admissible Prob optimal tc=" + std::to_string(tc) +
                              " mode=" + std::to_string(mode));
                if (ok) ++matched[mode];
            }
            freePath(p);
        }
    }

    std::cout << "Admissible Task 3 instances matched reference Dijkstra: "
              << "Mode1=" << matched[1] << "/" << admissibleTotal
              << ", Mode2=" << matched[2] << "/" << admissibleTotal << "\n";

    // Part B: General (potentially non-admissible) instances — record optimality gap
    int generalTotal = 60;
    int suboptCount[3] = {0, 0, 0};
    double maxGap[3] = {0.0, 0.0, 0.0};

    for (int tc = 0; tc < generalTotal; ++tc) {
        zeroMatrix(distM);
        zeroMatrix(riskM);
        zeroCoords(coords);
        int n = 15 + (tc % 20);
        double lambda = 10.0;

        for (int i = 0; i < n; ++i) {
            coords[i][0] = coordDist(rng);
            coords[i][1] = coordDist(rng);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j && probDist(rng) < 0.25) {
                    distM[i][j] = smallDDist(rng);
                    riskM[i][j] = riskDist(rng);
                }
            }
        }

        int s = static_cast<int>(rng() % n);
        int g = static_cast<int>(rng() % n);
        RefResult ref = refProbabilisticDijkstra(distM, riskM, s, g, lambda);

        for (int mode = 1; mode <= 2; ++mode) {
            PathNode* p = findProbabilisticPath(distM, riskM, coords, s, g, lambda, mode);
            if (!ref.found) {
                check(p == nullptr && std::fabs(g_lastPsafe) < EPS, "general Prob unreachable");
            } else {
                double totalG = -1.0, expPsafe = -1.0;
                bool valid = verifyProbPath(p, distM, riskM, coords, s, g, lambda, mode, totalG, expPsafe);
                check(valid && std::fabs(g_lastPsafe - expPsafe) < EPS, "general Prob valid & Psafe exact");
                double gap = totalG - ref.cost;
                check(gap >= -EPS, "general Prob cannot beat reference Dijkstra");
                if (gap > EPS) {
                    ++suboptCount[mode];
                    if (gap > maxGap[mode]) maxGap[mode] = gap;
                }
            }
            freePath(p);
        }
    }

    std::cout << std::fixed << std::setprecision(4)
              << "General (non-admissible) Task 3 gaps over " << generalTotal << " instances:\n"
              << "  Mode 1 (Manhattan):  " << suboptCount[1] << " suboptimal, maxGap=" << maxGap[1] << "\n"
              << "  Mode 2 (Euclidean):  " << suboptCount[2] << " suboptimal, maxGap=" << maxGap[2] << "\n";
}

} // namespace

int main() {
    testEdgeCasesAndLambdaTradeoff();
    testRandomAgainstReference();

    std::cout << "Task 3 Tests Summary: " << g_testsPassed << " passed, "
              << g_testsFailed << " failed.\n";
    return (g_testsFailed == 0) ? 0 : 1;
}
