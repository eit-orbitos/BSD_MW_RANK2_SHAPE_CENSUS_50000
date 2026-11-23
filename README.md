# Rank-two Mordell–Weil lattice shape census

**Study version V0.2 · Manuscript preprint v0.3.4**

This repository accompanies *A Reproducible Census of Rank-Two Mordell–Weil
Lattice Shapes for Elliptic Curves of Conductor Below 50,000* by **Toni
Mladenovski**, independent researcher, [ORCID 0009-0009-8343-662X](https://orcid.org/0009-0009-8343-662X).

It describes a fixed census from J. E. Cremona's `ecdata`: 23,612 rank-two curve
records and 18,886 predetermined representatives, one per isogeny class, for
conductors below 50,000. The mean reduced shape factor is approximately
0.953758660887 for all curves and 0.954750262395 for the representatives. The
normalized hyperbolic reference mean is 0.969197841658. These are descriptive
comparisons of a finite arithmetic population; no asymptotic distribution or
BSD proof is claimed.

## Start here

- [Manuscript PDF, version 0.3.4](manuscript/MW_Rank2_Census_v0_3_4_PREPRINT.pdf)
- [Reproduction instructions](REPRODUCING.md)
- [Stored summaries](recorded_results/summary_table.csv) and [selection of class representatives](recorded_results/population_B_selection.csv)
- [Citation instructions](CITATION.md) and [component licenses](licenses/LICENSE_SCOPE.md)

**Zenodo status: DOI_RESERVED_UNPUBLISHED.** The dataset DOI
`10.5281/zenodo.24484435` and preprint DOI `10.5281/zenodo.24484370` are reserved;
the corresponding Zenodo drafts are not yet published at this preparation step.
The v0.3.4 manuscript adds these identifiers and the GitHub repository link to
v0.3.3, with no scientific change. The prior v0.3.3 files remain as historical
copies. No human specialist peer review or journal acceptance is claimed. The
older BSD calculator is a separate work; its DOI does not identify this census.

## Contents and evidence

`code/` includes the ten scientific code files copied **byte for byte** from
V0.2, plus the standalone manuscript-figure script. The new
`scripts/fetch_inputs.py` helper retrieves the five exact upstream input blocks,
or obtains them from the canonical archive, and verifies their lengths and SHA-256
digests. It does not change the scientific computation.

The original V0.2 execution recorded 67 self-tests, 7 negative controls and
11 integration tests passing, with zero unresolved records. The AI-assisted
reproduction in a separate environment matched the compared scientific fields
across all 23,612 records. Both executions used PARI's height implementation;
they do not constitute verification by an independently implemented height
engine. Database rank and full-generator-basis assumptions are inherited.
See the manuscript's verification and AI-disclosure sections for the scope of
these checks. The files in `provenance/recorded_tests/` and the recorded exit
codes are retained historical outputs, not newly executed tests in this
repository-preparation step.

Historical review and reconciliation packages remain separate from this compact
repository. Their audits may depend on earlier manuscript filenames and review
fixtures; they are not represented here as standalone commands.

The large per-curve outputs and all five source blocks remain in the canonical
dataset archive, which is not duplicated in Git history:

```text
BSD_MW_SHAPE_STUDY_V0_2.zip
SHA-256 1338bdf8ddb50de5d9278a14e207ee3c3f2c9ab211d929b3556d7adf92a0d4fc
```

The source commit is `25cec5ecfec8b9f016eb1631ac633194c2bed39f` in
[JohnCremona/ecdata](https://github.com/JohnCremona/ecdata).
The exact original archive and this repository have distinct manifests.
Regenerated timestamps, execution times and complete output bytes need not
match those of the archived run.

## Licensing and attribution

| Component | License |
| --- | --- |
| Original scientific and support software | PolyForm Noncommercial 1.0.0 |
| Original manuscript, figures, documentation and author-held rights in the original result compilation | CC BY-NC 4.0 |
| Upstream ecdata inputs | Artistic License 2.0, unchanged |

The software is source available under the stated noncommercial terms.
These are separate component licenses; CC BY-NC is not an alternative software
license. Mathematical facts and third-party rights are not newly restricted or
relicensed by these notices. Preserve upstream attribution and notices.

**Funding:** This research received no external funding or grants.

**Competing interests:** The author declares no competing interests.

**AI assistance:** ChatGPT assisted with drafting and reconciliation; Claude
developed and executed the study pipeline and repairs; Kimi reran and reviewed
the outputs. These are AI-assisted technical checks, not human peer review.
The manuscript contains the full role and shared-engine disclosure. Toni
Mladenovski is the sole human author and retains responsibility for the work.
