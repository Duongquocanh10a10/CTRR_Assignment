#!/usr/bin/env bash
# Generates 5 standalone single-file C++ sources (one per task) under tools/coderunner/
# with NO local #include, ONLY allowed headers, and NO printPath, ready to paste into LMS CodeRunner.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

OUT_DIR="tools/coderunner"
mkdir -p "$OUT_DIR"

# Extract lines from submission/Algo.cpp for a 3-line banner section matching marker.
extract_section() {
    local marker="$1"
    awk -v m="$marker" '
        BEGIN { printing = 0 }
        { sub(/\r$/, "") }
        /^\/\/ ============================================================================/ {
            if (printing) { exit }
            banner = $0
            if (getline title > 0) {
                sub(/\r$/, "", title)
                if (getline closing > 0) {
                    sub(/\r$/, "", closing)
                }
                if (index(title, m) > 0) {
                    printing = 1
                    print banner
                    print title
                    print closing
                    next
                }
            }
            next
        }
        printing { print }
    ' submission/Algo.cpp
}

# Emit common preamble (allowed headers + PathNode struct + helpers + shared A* engine).
emit_common() {
    cat << 'EOF'
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

EOF
    extract_section "SHARED ENGINE"
}

# Task 1: findMeshPath
{
    emit_common
    extract_section "TASK 1 HELPERS"
    echo "} // namespace"
    echo ""
    extract_section "TASK 1:"
} > "$OUT_DIR/task1_findMeshPath.cpp"

# Task 2: findAUVPath
{
    emit_common
    extract_section "TASK 2 HELPERS"
    echo "} // namespace"
    echo ""
    extract_section "TASK 2:"
} > "$OUT_DIR/task2_findAUVPath.cpp"

# Task 3: findProbabilisticPath
{
    emit_common
    extract_section "TASK 3 HELPERS"
    echo "} // namespace"
    echo ""
    echo "double g_lastPsafe = 0.0;"
    echo ""
    extract_section "TASK 3:"
} > "$OUT_DIR/task3_findProbabilisticPath.cpp"

# Task 4: findGamePath
{
    emit_common
    extract_section "SHARED GRID CONVENTIONS"
    extract_section "TASK 4 HELPERS"
    echo "} // namespace"
    echo ""
    extract_section "TASK 4:"
} > "$OUT_DIR/task4_findGamePath.cpp"

# Task 5: findPCBPath
{
    emit_common
    extract_section "SHARED GRID CONVENTIONS"
    extract_section "TASK 4 HELPERS" | awk '/double heuristicGrid2D/,/^}/'
    echo ""
    extract_section "TASK 5 HELPERS" | grep -v 'g_lastPsafe' | grep -v '^// R10:'
    echo ""
    extract_section "TASK 5:"
} > "$OUT_DIR/task5_findPCBPath.cpp"

echo "Generated single-file CodeRunner sources in $OUT_DIR:"
ls -1 "$OUT_DIR"
