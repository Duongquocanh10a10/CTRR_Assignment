#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++11 -Wall -Wextra -Werror -O2}"
if [ "${1:-}" = "--asan" ]; then
    CXXFLAGS="-std=c++11 -Wall -Wextra -Werror -g -fsanitize=address,undefined"
fi

echo "=== 1. Checking allowed headers in submission/ ==="
ALLOWED="iostream|fstream|string|cmath|vector|algorithm|PathNode.h|Algo.h"
BAD_INCLUDES=$(grep -Hn '#include' submission/*.h submission/*.cpp | grep -Ev "#include[[:space:]]+[\"<]($ALLOWED)[\">]" || true)
if [ -n "$BAD_INCLUDES" ]; then
    echo "ERROR: Forbidden includes detected in submission/:"
    echo "$BAD_INCLUDES"
    exit 1
fi
echo "Header check PASSED (only allowed headers used)."

echo "=== 2. Checking printPath exclusion in PathNode.* and Algo.* ==="
if grep -Hn 'printPath' submission/PathNode.h submission/PathNode.cpp submission/Algo.h submission/Algo.cpp; then
    echo "ERROR: printPath must not appear in PathNode.* or Algo.*"
    exit 1
fi
echo "printPath exclusion PASSED."

echo "=== 3. Compiling and running submission/main.cpp ==="
$CXX $CXXFLAGS submission/PathNode.cpp submission/Algo.cpp submission/main.cpp -o tests/demo_main
./tests/demo_main > /dev/null
echo "submission/main.cpp compiled and ran cleanly."

echo "=== 4. Compiling and running all 5 task test suites ==="
for i in 1 2 3 4 5; do
    echo "--- Task $i ---"
    $CXX $CXXFLAGS submission/PathNode.cpp submission/Algo.cpp "tests/test_task${i}.cpp" -o "tests/test_task${i}"
    "./tests/test_task${i}"
done

echo "=== ALL CHECKS AND TESTS PASSED ==="
