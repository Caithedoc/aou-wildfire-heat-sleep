# Provenance

This file records what each document in `analysis_plan/` is, when it was written, and its SHA-256 checksum as
included in this repository. Checksums of the two 2026-10-12 addenda were also written into the project log
(`team_log.md`, kept with the authors) when the addenda were finalised, before the corresponding models were run.

| File | Date and role | SHA-256 of the file in this repository | Checksum recorded in the project log at finalisation |
|---|---|---|---|
| `analysis_plan/ANALYSIS_PLAN_20260826.md` | 2026-08-26 (written before data access; project file ANALYSIS_PLAN.md, last modified 2026-08-28 10:30) | `85a94b8bc37c7fa0361c1e12380abc4a053e6db93e0916beec3ae4032fa781ae` | not recorded |
| `analysis_plan/REVISION_PLAN_20260828.md` | 2026-08-28 (after the first full run; records the shift to the compound-first framing after the marginal smoke association proved non-robust) | `a1fc5238e22dd7fca5c2fa15139e4ab7168671f8af001e6f4c7267eb4401c57a` | not recorded |
| `analysis_plan/ANALYSIS_PLAN_ADDENDUM_SLEEPHEALTH_20261012.md` | 2026-10-12 (composite sleep-health outcomes, prespecified after coauthor review and before any composite model was fitted) | `160da447ec7dc7b366829ce3f144fd30856ef47b6bd5fe52676b7bc455c31652` | `160da447ec7dc7b366829ce3f144fd30856ef47b6bd5fe52676b7bc455c31652` |
| `analysis_plan/ANALYSIS_PLAN_ADDENDUM_SLEEPDOMAINS_20261012.md` | 2026-10-12 (sleep-domain secondary outcomes, prespecified after coauthor review and before any domain model was fitted) | `66fed30bbe6ca00dc48401dd372254a45fa269d581cf5b09ff51a690ac7157a8` | `92454e41...b674b5 (abbreviated in the project log)` |

Notes.

1. For `ANALYSIS_PLAN_ADDENDUM_SLEEPHEALTH_20261012.md` the checksum of the file included here matches the checksum
   recorded in the project log.
2. For `ANALYSIS_PLAN_ADDENDUM_SLEEPDOMAINS_20261012.md` the project log records an abbreviated checksum
   `92454e41...b674b5` at finalisation (before the domain models were launched on 2026-10-12). The file included here
   (last modified 2026-10-12 18:08 local time) has a different checksum, so it was edited on the same day after that
   record was made. The addendum's own rule permits later additions only to its post-hoc section (section 7); the
   pre-run copy was not retained separately, so the exact difference cannot be shown. We report this rather than
   claim byte-level immutability for that document.
3. The main plan (`ANALYSIS_PLAN_20260826.md`) and the revision plan were not checksummed at the time; their dates
   come from the documents themselves and from the project log.
4. The analysis-plan documents are in Chinese, as written. Section headings and the quantities they prespecify
   (exposure definitions, outcome definitions, model forms, sensitivity analyses) are summarised in English in
   `README.md`.

## Release edits

The scripts are the versions that produced the reported results, with one class of edit made for this release:
in the figure scripts, the hard-coded path to the authors' local project folder (`BASE = "/Users/.../CZ_AoU_WildfireSleep_08_26_26_V1"`)
was replaced by the repository root (`BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))`),
with `import os` added where it was missing. No analytic code was changed. Scripts that belong to other projects
(a feasibility-count script for a separate study) and archived earlier figure drafts were not included.

## File checksums

`MANIFEST.sha256` lists the SHA-256 checksum of every file in the repository at release.
