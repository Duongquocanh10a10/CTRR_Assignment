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

void zeroGrid(int g[100][100]) {
    std::memset(g, 0, sizeof(int) * 100 * 100);
}

bool dirLabelToDelta(const std::string& label, int& dx, int& dy) {
    if (label == "Up") { dx = -1; dy = 0; return true; }
    if (label == "Down") { dx = 1; dy = 0; return true; }
    if (label == "Left") { dx = 0; dy = -1; return true; }
    if (label == "Right") { dx = 0; dy = 1; return true; }
    if (label == "Up-Left") { dx = -1; dy = -1; return true; }
    if (label == "Up-Right") { dx = -1; dy = 1; return true; }
    if (label == "Down-Left") { dx = 1; dy = -1; return true; }
    if (label == "Down-Right") { dx = 1; dy = 1; return true; }
    return false;
}

double expectedH(int x, int y, int gx, int gy, int mode) {
    double dx = std::fabs(static_cast<double>(x - gx));
    double dy = std::fabs(static_cast<double>(y - gy));
    return (mode == 2) ? std::max(dx, dy) : (dx + dy);
}

bool verifyGamePath(PathNode* head,
                    int mapGrid[100][100], int m, int n,
                    int startX, int startY, int goalX, int goalY, int mode,
                    double& outTotalG) {
    if (head == nullptr) return false;
    int effMode = (mode == 1 || mode == 2) ? mode : 1;

    if (head->name != "Start") return false;
    if (std::fabs(head->g) > EPS) return false;

    int x = startX, y = startY;
    double accumG = 0.0;

    PathNode* cur = head;
    bool first = true;
    while (cur != nullptr) {
        if (!first) {
            int dx = 0, dy = 0;
            if (!dirLabelToDelta(cur->name, dx, dy)) return false;
            int nx = x + dx;
            int ny = y + dy;
            if (nx < 0 || nx >= m || ny < 0 || ny >= n) return false;
            if (mapGrid[nx][ny] == 1) return false;
            bool diag = (dx != 0 && dy != 0);
            if (diag) {
                if (mapGrid[x + dx][y] == 1 || mapGrid[x][y + dy] == 1) return false;
            }
            double base = (mapGrid[nx][ny] == 2) ? 2.0 : 1.0;
            accumG += diag ? base * 1.5 : base;
            x = nx;
            y = ny;
        }
        double hExp = expectedH(x, y, goalX, goalY, effMode);
        if (std::fabs(cur->g - accumG) > EPS) return false;
        if (std::fabs(cur->h - hExp) > EPS) return false;
        if (std::fabs(cur->f - (cur->g + cur->h)) > EPS) return false;

        first = false;
        cur = cur->next;
    }

    if (x != goalX || y != goalY) return false;
    outTotalG = accumG;
    return true;
}

void testEdgeCases() {
    int grid[100][100];
    zeroGrid(grid);

    // 1. start == goal (flat & swamp)
    {
        grid[2][2] = 2; // swamp
        PathNode* p = findGamePath(grid, 5, 5, 2, 2, 2, 2, 2);
        check(p != nullptr, "start == goal on swamp returns non-null");
        if (p != nullptr) {
            check(p->name == "Start" && std::fabs(p->g) < EPS &&
                      std::fabs(p->h) < EPS && std::fabs(p->f) < EPS && p->next == nullptr,
                  "start == goal single node Start with f=g=h=0");
            freePath(p);
        }
    }

    // 2. Start or goal on obstacle, out of bounds
    {
        zeroGrid(grid);
        grid[0][0] = 1;
        check(findGamePath(grid, 5, 5, 0, 0, 4, 4, 2) == nullptr, "start on obstacle -> nullptr");
        grid[0][0] = 0;
        grid[4][4] = 1;
        check(findGamePath(grid, 5, 5, 0, 0, 4, 4, 2) == nullptr, "goal on obstacle -> nullptr");
        check(findGamePath(grid, 0, 5, 0, 0, 0, 1, 2) == nullptr, "m=0 -> nullptr");
        check(findGamePath(grid, 5, 5, -1, 0, 2, 2, 2) == nullptr, "startX<0 -> nullptr");
    }

    // 3. Corner-cutting blocked:
    // (0,1)=1 and (1,0)=1 block diagonal from (0,0) to (1,1) -> no path at all!
    {
        zeroGrid(grid);
        grid[0][1] = 1;
        grid[1][0] = 1;
        PathNode* p = findGamePath(grid, 2, 2, 0, 0, 1, 1, 2);
        check(p == nullptr, "corner-cutting blocked completely prevents (0,0)->(1,1)");
    }

    // 4. Single corner obstacle blocks diagonal, forcing 2 cardinal steps (cost 2.0 instead of 1.5)
    {
        zeroGrid(grid);
        grid[0][1] = 1;
        PathNode* p = findGamePath(grid, 2, 2, 0, 0, 1, 1, 2);
        double totalG = -1.0;
        check(verifyGamePath(p, grid, 2, 2, 0, 0, 1, 1, 2, totalG), "single corner obstacle path valid");
        check(std::fabs(totalG - 2.0) < EPS, "forced around corner costs 2.0 (Down then Right)");
        freePath(p);
    }

    // 5. Open diagonal step costs 1.5 with label "Down-Right"
    {
        zeroGrid(grid);
        PathNode* p = findGamePath(grid, 2, 2, 0, 0, 1, 1, 2);
        double totalG = -1.0;
        check(verifyGamePath(p, grid, 2, 2, 0, 0, 1, 1, 2, totalG), "open diagonal valid");
        check(std::fabs(totalG - 1.5) < EPS, "open diagonal costs 1.5");
        if (p != nullptr && p->next != nullptr) {
            check(p->next->name == "Down-Right", "open diagonal label is Down-Right");
        }
        freePath(p);
    }

    // 6. Swamp vs detour trade-off:
    // Grid 3x3: start (1,0), goal (1,2).
    // Direct through (1,1)=swamp: (1,0)->(1,1)->(1,2) costs 2.0 + 1.0 = 3.0.
    // Diagonal detour via (0,1)=flat: (1,0)->(0,1) [Up-Right, 1.5] -> (1,2) [Down-Right, 1.5] = 3.0.
    // Now make (1,1) swamp and (1,2) flat, and check 3x4 with two swamps (1,1)=2, (1,2)=2:
    // Direct through 2 swamps: 2 + 2 + 1 = 5.0.
    // Detour via row 0 flat: (1,0)->(0,1)[1.5]->(0,2)[1.0]->(1,3)[1.5] = 4.0!
    {
        zeroGrid(grid);
        grid[1][1] = 2;
        grid[1][2] = 2;
        PathNode* p = findGamePath(grid, 3, 4, 1, 0, 1, 3, 2);
        double totalG = -1.0;
        check(verifyGamePath(p, grid, 3, 4, 1, 0, 1, 3, 2, totalG), "swamp detour valid");
        check(std::fabs(totalG - 4.0) < EPS, "swamp detour avoids 2 swamps (4.0 < 5.0)");
        freePath(p);
    }

    // 7. Max 100x100 grid test
    {
        zeroGrid(grid);
        for (int i = 0; i < 95; ++i) {
            grid[i][50] = 1; // tall wall forcing a detour
        }
        RefResult ref = refGameDijkstra(grid, 100, 100, 0, 0, 0, 99);
        PathNode* p = findGamePath(grid, 100, 100, 0, 0, 0, 99, 2);
        double totalG = -1.0;
        check(verifyGamePath(p, grid, 100, 100, 0, 0, 0, 99, 2, totalG), "100x100 grid valid");
        check(ref.found && std::fabs(totalG - ref.cost) < EPS, "100x100 grid optimal in Mode 2");
        freePath(p);
    }
}

void testRandomAgainstReference() {
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> uDist(0.0, 1.0);

    int grid[100][100];
    int totalInstances = 100;
    int mode2Matched = 0;
    int mode1Subopt = 0;
    double mode1MaxGap = 0.0;

    for (int tc = 0; tc < totalInstances; ++tc) {
        zeroGrid(grid);
        int m = 6 + (tc % 15);
        int n = 6 + ((tc * 3) % 15);
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                double r = uDist(rng);
                if (r < 0.22) grid[i][j] = 1;       // obstacle
                else if (r < 0.45) grid[i][j] = 2;  // swamp
                else grid[i][j] = 0;                // flat
            }
        }
        int sx = static_cast<int>(rng() % m);
        int sy = static_cast<int>(rng() % n);
        int gx = static_cast<int>(rng() % m);
        int gy = static_cast<int>(rng() % n);
        grid[sx][sy] = 0;
        grid[gx][gy] = 0;

        RefResult ref = refGameDijkstra(grid, m, n, sx, sy, gx, gy);

        // Mode 2 (Chebyshev): admissible -> must match reference Dijkstra exactly
        {
            PathNode* p = findGamePath(grid, m, n, sx, sy, gx, gy, 2);
            if (!ref.found) {
                check(p == nullptr, "random Game Mode 2 unreachable");
                if (p == nullptr) ++mode2Matched;
            } else {
                double totalG = -1.0;
                bool ok = verifyGamePath(p, grid, m, n, sx, sy, gx, gy, 2, totalG) &&
                          (std::fabs(totalG - ref.cost) < EPS);
                check(ok, "random Game Mode 2 optimal tc=" + std::to_string(tc));
                if (ok) ++mode2Matched;
            }
            freePath(p);
        }

        // Mode 1 (Manhattan): non-admissible with diagonal cost 1.5 -> verify validity & record gap
        {
            PathNode* p = findGamePath(grid, m, n, sx, sy, gx, gy, 1);
            if (!ref.found) {
                check(p == nullptr, "random Game Mode 1 unreachable");
            } else {
                double totalG = -1.0;
                bool valid = verifyGamePath(p, grid, m, n, sx, sy, gx, gy, 1, totalG);
                check(valid, "random Game Mode 1 valid tc=" + std::to_string(tc));
                double gap = totalG - ref.cost;
                check(gap >= -EPS, "random Game Mode 1 cannot beat reference Dijkstra");
                if (gap > EPS) {
                    ++mode1Subopt;
                    if (gap > mode1MaxGap) mode1MaxGap = gap;
                }
            }
            freePath(p);
        }
    }

    std::cout << "Random Task 4 Mode 2 (Chebyshev, admissible) matched reference Dijkstra: "
              << mode2Matched << "/" << totalInstances << "\n";
    std::cout << std::fixed << std::setprecision(4)
              << "Random Task 4 Mode 1 (Manhattan, non-admissible) gaps: "
              << mode1Subopt << "/" << totalInstances << " suboptimal, maxGap=" << mode1MaxGap << "\n";
}

} // namespace

int main() {
    testEdgeCases();
    testRandomAgainstReference();

    std::cout << "Task 4 Tests Summary: " << g_testsPassed << " passed, "
              << g_testsFailed << " failed.\n";
    return (g_testsFailed == 0) ? 0 : 1;
}
