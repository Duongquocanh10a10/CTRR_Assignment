#ifndef REF_DIJKSTRA_H
#define REF_DIJKSTRA_H

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

// Independent brute-force Dijkstra solver (no heuristic) used as a ground-truth
// reference oracle for verifying optimal path costs across all 5 tasks.
struct RefResult {
    bool found;
    double cost;
    std::vector<int> states;
};

inline RefResult refDijkstra(
    int numStates,
    int startState,
    const std::function<bool(int)>& isGoal,
    const std::function<void(int, const std::function<void(int, double)>&)>& forEachNeighbor) {
    const double INF = 1e18;
    const double EPS = 1e-12;

    if (numStates <= 0 || startState < 0 || startState >= numStates) {
        return {false, INF, {}};
    }

    std::vector<double> dist(numStates, INF);
    std::vector<int> parent(numStates, -1);

    using Entry = std::pair<double, int>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> pq;

    dist[startState] = 0.0;
    pq.push({0.0, startState});

    int goalState = -1;
    while (!pq.empty()) {
        double d = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (d > dist[u] + EPS) {
            continue;
        }
        if (isGoal(u)) {
            goalState = u;
            break;
        }

        forEachNeighbor(u, [&](int v, double w) {
            double nd = d + w;
            if (nd < dist[v] - EPS) {
                dist[v] = nd;
                parent[v] = u;
                pq.push({nd, v});
            }
        });
    }

    if (goalState == -1) {
        return {false, INF, {}};
    }

    std::vector<int> path;
    for (int s = goalState; s != -1; s = parent[s]) {
        path.push_back(s);
    }
    std::reverse(path.begin(), path.end());
    return {true, dist[goalState], path};
}

// Reference solver specifically for Task 1 (Mesh Network with Congestion Penalty).
inline RefResult refMeshDijkstra(double weightMatrix[100][100], int startNode, int goalNode) {
    if (startNode < 0 || startNode >= 100 || goalNode < 0 || goalNode >= 100) {
        return {false, 1e18, {}};
    }

    int deg[100] = {0};
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            if (i != j && weightMatrix[i][j] > 0.0) {
                ++deg[i];
            }
        }
    }

    return refDijkstra(
        100,
        startNode,
        [&](int s) { return s == goalNode; },
        [&](int u, const std::function<void(int, double)>& visit) {
            for (int v = 0; v < 100; ++v) {
                if (v != u && weightMatrix[u][v] > 0.0) {
                    visit(v, weightMatrix[u][v] + 0.1 * static_cast<double>(deg[v]));
                }
            }
        });
}

// Reference solver specifically for Task 2 (AUV 3D Navigation with Depth Penalty).
inline RefResult refAUVDijkstra(double weightMatrix[100][100],
                                int coords[100][3],
                                int startPoint,
                                int goalPoint) {
    if (startPoint < 0 || startPoint >= 100 || goalPoint < 0 || goalPoint >= 100) {
        return {false, 1e18, {}};
    }

    return refDijkstra(
        100,
        startPoint,
        [&](int s) { return s == goalPoint; },
        [&](int u, const std::function<void(int, double)>& visit) {
            for (int v = 0; v < 100; ++v) {
                if (v != u && weightMatrix[u][v] > 0.0) {
                    double dz = std::fabs(static_cast<double>(coords[v][2] - coords[u][2]));
                    visit(v, weightMatrix[u][v] * (1.0 + dz / 100.0));
                }
            }
        });
}

// Reference solver specifically for Task 3 (Probabilistic Pathfinding under Hazard Risk).
inline RefResult refProbabilisticDijkstra(double distMatrix[100][100],
                                          double riskMatrix[100][100],
                                          int startPoint,
                                          int goalPoint,
                                          double lambda) {
    if (startPoint < 0 || startPoint >= 100 || goalPoint < 0 || goalPoint >= 100) {
        return {false, 1e18, {}};
    }

    return refDijkstra(
        100,
        startPoint,
        [&](int s) { return s == goalPoint; },
        [&](int u, const std::function<void(int, double)>& visit) {
            for (int v = 0; v < 100; ++v) {
                if (v != u && distMatrix[u][v] > 0.0) {
                    visit(v, distMatrix[u][v] + lambda * riskMatrix[u][v]);
                }
            }
        });
}

// Reference solver specifically for Task 4 (8-dir Multi-Terrain Game Grid).
inline RefResult refGameDijkstra(int mapGrid[100][100], int m, int n,
                                 int startX, int startY, int goalX, int goalY) {
    if (m <= 0 || m > 100 || n <= 0 || n > 100) return {false, 1e18, {}};
    if (startX < 0 || startX >= m || startY < 0 || startY >= n) return {false, 1e18, {}};
    if (goalX < 0 || goalX >= m || goalY < 0 || goalY >= n) return {false, 1e18, {}};
    if (mapGrid[startX][startY] == 1 || mapGrid[goalX][goalY] == 1) return {false, 1e18, {}};

    const int dx[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
    const int dy[8] = {0, 0, -1, 1, -1, 1, -1, 1};
    int goalCell = goalX * n + goalY;

    return refDijkstra(
        m * n,
        startX * n + startY,
        [&](int s) { return s == goalCell; },
        [&](int u, const std::function<void(int, double)>& visit) {
            int x = u / n;
            int y = u % n;
            for (int d = 0; d < 8; ++d) {
                int nx = x + dx[d];
                int ny = y + dy[d];
                if (nx < 0 || nx >= m || ny < 0 || ny >= n) continue;
                if (mapGrid[nx][ny] == 1) continue;
                bool diag = (dx[d] != 0 && dy[d] != 0);
                if (diag) {
                    if (mapGrid[x + dx[d]][y] == 1 || mapGrid[x][y + dy[d]] == 1) continue;
                }
                double base = (mapGrid[nx][ny] == 2) ? 2.0 : 1.0;
                visit(nx * n + ny, diag ? base * 1.5 : base);
            }
        });
}

// Reference solver specifically for Task 5 (4-dir PCB Wire Routing with Bending Penalty).
inline RefResult refPCBDijkstra(int pcbGrid[100][100], int m, int n,
                                int startX, int startY, int exitX, int exitY) {
    if (m <= 0 || n <= 0 || m * n > 100) return {false, 1e18, {}};
    if (startX < 0 || startX >= m || startY < 0 || startY >= n) return {false, 1e18, {}};
    if (exitX < 0 || exitX >= m || exitY < 0 || exitY >= n) return {false, 1e18, {}};
    if (pcbGrid[startX][startY] != 0 || pcbGrid[exitX][exitY] != 0) return {false, 1e18, {}};

    const int dx[4] = {-1, 1, 0, 0};
    const int dy[4] = {0, 0, -1, 1};
    int exitCell = exitX * n + exitY;

    return refDijkstra(
        m * n * 5,
        (startX * n + startY) * 5 + 4,
        [&](int s) { return (s / 5) == exitCell; },
        [&](int s, const std::function<void(int, double)>& visit) {
            int cell = s / 5;
            int prevDir = s % 5;
            int x = cell / n;
            int y = cell % n;
            for (int d = 0; d < 4; ++d) {
                int nx = x + dx[d];
                int ny = y + dy[d];
                if (nx < 0 || nx >= m || ny < 0 || ny >= n) continue;
                if (pcbGrid[nx][ny] != 0) continue;
                double bend = (prevDir != 4 && prevDir != d) ? 0.5 : 0.0;
                visit((nx * n + ny) * 5 + d, 1.0 + bend);
            }
        });
}

#endif // REF_DIJKSTRA_H
