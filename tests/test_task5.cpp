#include "../submission/Algo.h"
#include "../submission/PathNode.h"
#include "reference/ref_dijkstra.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
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

void zeroGrid(int g[100][100]) {
    std::memset(g, 0, sizeof(int) * 100 * 100);
}

void dirtyMatrix(double m[100][100]) {
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            m[i][j] = 999.0;
        }
    }
}

bool parseCoord(const std::string& name, int& x, int& y) {
    return std::sscanf(name.c_str(), "(%d,%d)", &x, &y) == 2;
}

double expectedH(int x, int y, int ex, int ey, int mode) {
    double dx = std::fabs(static_cast<double>(x - ex));
    double dy = std::fabs(static_cast<double>(y - ey));
    return (mode == 2) ? std::max(dx, dy) : (dx + dy);
}

bool verifyWeightMatrix(int pcbGrid[100][100], int m, int n, double w[100][100]) {
    int numCells = m * n;
    for (int u = 0; u < numCells; ++u) {
        int ux = u / n, uy = u % n;
        for (int v = 0; v < numCells; ++v) {
            int vx = v / n, vy = v % n;
            bool adj = (std::abs(ux - vx) + std::abs(uy - vy) == 1);
            bool expectedEdge = adj && (pcbGrid[ux][uy] == 0) && (pcbGrid[vx][vy] == 0);
            double expW = expectedEdge ? 1.0 : 0.0;
            if (std::fabs(w[u][v] - expW) > EPS) {
                return false;
            }
        }
    }
    return true;
}

bool verifyPCBPath(PathNode* head,
                   int pcbGrid[100][100], int m, int n,
                   int startX, int startY, int exitX, int exitY, int mode,
                   double& outTotalG) {
    if (head == nullptr) return false;
    int effMode = (mode == 1 || mode == 2) ? mode : 1;

    std::vector<PathNode*> nodes;
    for (PathNode* cur = head; cur != nullptr; cur = cur->next) {
        nodes.push_back(cur);
    }

    int px = -1, py = -1;
    int prevDir = -1; // -1 = none
    double accumG = 0.0;

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        int x = -1, y = -1;
        if (!parseCoord(nodes[i]->name, x, y)) return false;
        if (x < 0 || x >= m || y < 0 || y >= n) return false;
        if (pcbGrid[x][y] != 0) return false;

        if (i == 0) {
            if (x != startX || y != startY) return false;
        } else {
            int dx = x - px;
            int dy = y - py;
            if (std::abs(dx) + std::abs(dy) != 1) return false;
            int dir = (dx == -1) ? 0 : (dx == 1) ? 1 : (dy == -1) ? 2 : 3;
            double bend = (prevDir != -1 && prevDir != dir) ? 0.5 : 0.0;
            accumG += 1.0 + bend;
            prevDir = dir;
        }

        double hExp = expectedH(x, y, exitX, exitY, effMode);
        if (std::fabs(nodes[i]->g - accumG) > EPS) return false;
        if (std::fabs(nodes[i]->h - hExp) > EPS) return false;
        if (std::fabs(nodes[i]->f - (nodes[i]->g + nodes[i]->h)) > EPS) return false;

        px = x;
        py = y;
    }

    if (px != exitX || py != exitY) return false;
    outTotalG = accumG;
    return true;
}

void testEdgeCasesAndBendingTradeoff() {
    int grid[100][100];
    double w[100][100];

    // 1. start == exit -> single node "(x,y)" with f=g=h=0, weightMatrix populated
    {
        zeroGrid(grid);
        dirtyMatrix(w);
        PathNode* p = findPCBPath(grid, 4, 4, 1, 2, 1, 2, w, 1);
        check(verifyWeightMatrix(grid, 4, 4, w), "weightMatrix populated on start==exit");
        check(p != nullptr, "start == exit returns non-null");
        if (p != nullptr) {
            check(p->name == "(1,2)" && std::fabs(p->g) < EPS &&
                      std::fabs(p->h) < EPS && std::fabs(p->f) < EPS && p->next == nullptr,
                  "start == exit single node (1,2) with f=g=h=0");
            freePath(p);
        }
    }

    // 2. m * n > 100 -> nullptr (A6)
    {
        zeroGrid(grid);
        dirtyMatrix(w);
        check(findPCBPath(grid, 11, 10, 0, 0, 5, 5, w, 1) == nullptr, "m*n=110 > 100 -> nullptr");
        check(findPCBPath(grid, 0, 5, 0, 0, 0, 1, w, 1) == nullptr, "m=0 -> nullptr");
    }

    // 3. Start or exit on obstacle -> nullptr (weightMatrix still populated)
    {
        zeroGrid(grid);
        dirtyMatrix(w);
        grid[0][0] = 1;
        check(findPCBPath(grid, 4, 4, 0, 0, 3, 3, w, 1) == nullptr, "start on obstacle -> nullptr");
        check(verifyWeightMatrix(grid, 4, 4, w), "weightMatrix still populated when start blocked");
    }

    // 4. Disconnected / walled off -> nullptr
    {
        zeroGrid(grid);
        dirtyMatrix(w);
        for (int j = 0; j < 5; ++j) grid[2][j] = 1;
        check(findPCBPath(grid, 5, 5, 0, 0, 4, 4, w, 1) == nullptr, "walled off -> nullptr");
    }

    // 5. Bending makes a LONGER straight route win!
    // 6x5 grid, Start (0,0), Exit (4,4).
    // Corridor A (tight staircase): (0,0)->(0,1)->(1,1)->(1,2)->(2,2)->(2,3)->(3,3)->(3,4)->(4,4)
    //   8 steps, 7 bends => cost = 8 + 7 * 0.5 = 11.5.
    // Corridor B (straight outer route): (0,0)->(5,0) [5 Down] -> (5,4) [4 Right] -> (4,4) [1 Up]
    //   10 steps, 2 bends => cost = 10 + 2 * 0.5 = 11.0!
    // Block everything else so only Corridor A and Corridor B exist.
    {
        zeroGrid(grid);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 5; ++j) {
                grid[i][j] = 1; // start all blocked
            }
        }
        // Carve Corridor A (staircase)
        int corrA[][2] = {{0,0},{0,1},{1,1},{1,2},{2,2},{2,3},{3,3},{3,4},{4,4}};
        for (auto& cell : corrA) grid[cell[0]][cell[1]] = 0;
        // Carve Corridor B (outer straight)
        for (int r = 0; r <= 5; ++r) grid[r][0] = 0;
        for (int c = 0; c <= 4; ++c) grid[5][c] = 0;
        grid[4][4] = 0;

        dirtyMatrix(w);
        for (int mode = 1; mode <= 2; ++mode) {
            PathNode* p = findPCBPath(grid, 6, 5, 0, 0, 4, 4, w, mode);
            double totalG = -1.0;
            check(verifyWeightMatrix(grid, 6, 5, w), "bending tradeoff weightMatrix valid");
            check(verifyPCBPath(p, grid, 6, 5, 0, 0, 4, 4, mode, totalG),
                  "bending tradeoff path valid mode " + std::to_string(mode));
            check(std::fabs(totalG - 11.0) < EPS,
                  "longer straight route (10 steps, 2 bends = 11.0) beats staircase (11.5)");
            // Verify path actually took (1,0) instead of (0,1)
            if (p != nullptr && p->next != nullptr) {
                check(p->next->name == "(1,0)", "bending tradeoff first move goes Down to (1,0)");
            }
            freePath(p);
        }
    }

    // 6. Max 10x10 (100 cells) grid
    {
        zeroGrid(grid);
        dirtyMatrix(w);
        for (int r = 0; r < 8; ++r) grid[r][4] = 1;
        RefResult ref = refPCBDijkstra(grid, 10, 10, 0, 0, 0, 9);
        for (int mode = 1; mode <= 2; ++mode) {
            PathNode* p = findPCBPath(grid, 10, 10, 0, 0, 0, 9, w, mode);
            double totalG = -1.0;
            check(verifyWeightMatrix(grid, 10, 10, w), "10x10 weightMatrix valid");
            check(verifyPCBPath(p, grid, 10, 10, 0, 0, 0, 9, mode, totalG),
                  "10x10 PCB path valid mode " + std::to_string(mode));
            check(ref.found && std::fabs(totalG - ref.cost) < EPS,
                  "10x10 PCB optimal mode " + std::to_string(mode));
            freePath(p);
        }
    }
}

void testRandomAgainstReference() {
    std::mt19937 rng(9001);
    std::uniform_real_distribution<double> uDist(0.0, 1.0);

    int grid[100][100];
    double w[100][100];
    int totalInstances = 100;
    int matched[3] = {0, 0, 0};

    for (int tc = 0; tc < totalInstances; ++tc) {
        zeroGrid(grid);
        int m = 4 + (tc % 7); // 4..10
        int n = 4 + ((tc * 3) % 7); // 4..10 => m*n <= 100
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                grid[i][j] = (uDist(rng) < 0.25) ? 1 : 0;
            }
        }
        int sx = static_cast<int>(rng() % m);
        int sy = static_cast<int>(rng() % n);
        int ex = static_cast<int>(rng() % m);
        int ey = static_cast<int>(rng() % n);
        grid[sx][sy] = 0;
        grid[ex][ey] = 0;

        RefResult ref = refPCBDijkstra(grid, m, n, sx, sy, ex, ey);

        for (int mode = 1; mode <= 2; ++mode) {
            dirtyMatrix(w);
            PathNode* p = findPCBPath(grid, m, n, sx, sy, ex, ey, w, mode);
            check(verifyWeightMatrix(grid, m, n, w), "random PCB weightMatrix valid");
            if (!ref.found) {
                check(p == nullptr, "random PCB unreachable mode " + std::to_string(mode));
                if (p == nullptr) ++matched[mode];
            } else {
                double totalG = -1.0;
                bool ok = verifyPCBPath(p, grid, m, n, sx, sy, ex, ey, mode, totalG) &&
                          (std::fabs(totalG - ref.cost) < EPS);
                check(ok, "random PCB optimal tc=" + std::to_string(tc) +
                              " mode=" + std::to_string(mode));
                if (ok) ++matched[mode];
            }
            freePath(p);
        }
    }

    std::cout << "Random Task 5 instances matched reference Dijkstra: "
              << "Mode1=" << matched[1] << "/" << totalInstances
              << ", Mode2=" << matched[2] << "/" << totalInstances << "\n";
}

} // namespace

int main() {
    testEdgeCasesAndBendingTradeoff();
    testRandomAgainstReference();

    std::cout << "Task 5 Tests Summary: " << g_testsPassed << " passed, "
              << g_testsFailed << " failed.\n";
    return (g_testsFailed == 0) ? 0 : 1;
}
