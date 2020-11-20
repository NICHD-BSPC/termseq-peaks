#!/bin/bash
set -e
termseq_peaks reverse_*.gz --peaks out.bed
diff out.bed expected_output.bed
