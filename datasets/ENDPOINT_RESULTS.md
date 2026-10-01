# T4.2 held-out endpoint results

Generated from `datasets/heldout/endpoint_results/*/endpoints.json` by
`verification/validate_endpoints.py`, against the spec in `datasets/ENDPOINTS.md`
(signed 2026-10-01, amended 6.1 A1-A6). Raw per-endpoint JSON lives under
`datasets/heldout/`, which is gitignored with the data it describes.

| stratum | run | endpoint | verdict | key observed |
|---|---|---|---|---|
| cross_lab | ERR10015758 | E1 | PASS | dominant 0.4005 |
| cross_lab | ERR10015758 | E3 | FAIL | r=0.6269 |
| cross_lab | ERR10015758 | E5 | PASS | known 0.7063 |
| cross_lab | ERR10015758 | E2 | NOT_EVALUATED |  |
| cross_lab | ERR10015758 | E4 | NOT_EVALUATED |  |
| cross_lab | ERR10015758 | E6 | NOT_EVALUATED |  |
| cross_chemistry | ERR10229623 | E1 | INCONCLUSIVE | unassigned 0.1546 |
| cross_chemistry | ERR10229623 | E3 | PASS | r=0.9475 |
| cross_chemistry | ERR10229623 | E5 | FAIL | known 0.4109 |
| cross_chemistry | ERR10229623 | E2 | NOT_EVALUATED |  |
| cross_chemistry | ERR10229623 | E4 | NOT_EVALUATED |  |
| cross_chemistry | ERR10229623 | E6 | NOT_EVALUATED |  |
| cross_organism | SRR1177982 | E1 | INCONCLUSIVE | unassigned 0.0009 |
| cross_organism | SRR1177982 | E3 | FAIL | r=0.67 |
| cross_organism | SRR1177982 | E5 | INCONCLUSIVE | known 0.0257 |
| cross_organism | SRR1177982 | E2 | NOT_EVALUATED |  |
| cross_organism | SRR1177982 | E4 | NOT_EVALUATED |  |
| cross_organism | SRR1177982 | E6 | NOT_EVALUATED |  |

## Reading these

E3 and E5 are NOT tests of the port. `infer_experiment` on the cross-lab BAM is
byte-identical between the port and upstream RSeQC (0.2024 / 0.3971 / 0.4005),
so a divergence from the dev reference is a property of the library, not a defect.

- **E3** passes on BGISEQ-500 (r=0.948) and fails on HiSeq 2500 (0.627) and rat
  (0.670): gene-body skew shape tracks library prep, so one threshold across
  chemistries measures the chemistry.
- **E5** passes on the Charité run (0.706) and fails on the BGI run (0.411) against
  the SAME GENCODE annotation on the SAME contigs, so 0.411 is a real failure and
  not an annotation artifact. Rat is INCONCLUSIVE for a different reason: its
  annotation covers 32% of the index, against 503% for human GENCODE.
- **E1** is INCONCLUSIVE wherever a contig-subset index confounds the unassigned
  fraction, and evaluates as a strand-BALANCE check where ENA metadata says
  `library_selection=PCR`.
- **E2, E4, E6** are not evaluable on any held-out stratum: no ERCC spike-ins
  exist in any of the three runs (E4), the degradation re-alignment was not run
  (E2), and no single-cell data exists (E6).
