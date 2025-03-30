#!/bin/bash

echo "🔍 Running flake8..."
flake8 clustermolepy/
FLAKE8_EXIT=$?

echo "📦 Checking import order with isort..."
isort clustermolepy/ --check
ISORT_EXIT=$?

echo "🖤 Checking code formatting with black..."
black clustermolepy/ --check
BLACK_EXIT=$?

# Summary
echo ""
if [[ $FLAKE8_EXIT -eq 0 && $ISORT_EXIT -eq 0 && $BLACK_EXIT -eq 0 ]]; then
    echo "✅ All checks passed!"
    exit 0
else
    echo "❌ Some checks failed:"
    [[ $FLAKE8_EXIT -ne 0 ]] && echo " - flake8"
    [[ $ISORT_EXIT -ne 0 ]] && echo " - isort"
    [[ $BLACK_EXIT -ne 0 ]] && echo " - black"
    exit 1
fi

