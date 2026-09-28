#!/bin/bash
# mock htseq-count: ignores its arguments
cat "$(dirname "$0")/htseq_counts.txt"
