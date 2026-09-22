# Restart runbook -- shutting down and coming back with context

**2026-07-28.** Written for a machine restart during Block B (figure phase). Applies to any
future one.

---

## The short version

**No analysis script needs re-running.** `results/*.rds` are files on disk, not session
state, and a reboot cannot touch them. Every input the whole `figures/` layer reads is
present and newer than the script that wrote it (checked 2026-07-28). What a restart costs
is the R session's memory and the terminal; what it *risks* is Google Drive.

Coming back is two commands:

```r
source(here::here("scripts", "99_session_bootstrap.R"))   # in Positron
```
```bash
./paper/clean_icon_files.sh repo                          # in a terminal, FIRST
```

---

## Before shutdown

1. **Confirm clean and pushed.**
   ```bash
   git status --short && git status -sb | head -1
   ```
   The only expected untracked file is `docs/library_reference/gray_chea_mito_tf_shortlist.csv`
   (a harmless duplicate of the tracked copy under `data/genesets_from_library/`; leave it,
   and never `git add -A docs`).

2. **Let Google Drive finish syncing.** Wait for the menu-bar item to say everything is up
   to date before powering off. A restart in the middle of a sync is how the tree ends up in
   a strange state.

3. **Optional, and cheap insurance.** `results/` is gitignored, so it exists in exactly one
   place -- and two of its objects cannot be regenerated from the trunk (see hazards below).
   A snapshot on local disk, off the Drive mount:
   ```bash
   tar czf ~/myc_results_$(date +%F).tgz -C /Users/gs/G/data/MK_myc_2022/myc_mouse results
   ```

4. **Nothing to save from the R session.** Do not save the workspace. Everything the figure
   layer needs is already in `results/`; a saved `.RData` is a liability, not a backup (see
   hazards).

---

## After restart

1. **Let Drive mount and settle** before touching the repo. It re-syncs the whole tree on
   boot, and that is exactly when it stamps `Icon\r` stubs everywhere.

2. **Sweep the stubs, before the first git command.**
   ```bash
   cd /Users/gs/G/data/MK_myc_2022/myc_mouse
   ./paper/clean_icon_files.sh repo
   ```
   This is the whole repo including `.git`. It deletes only files whose name is exactly
   `Icon` + carriage return **and** which are zero bytes -- no legitimate git file can be
   both, so it cannot damage the repository.

3. **Verify git.**
   ```bash
   git rev-parse HEAD && git status --short && git fetch origin
   ```
   If `fetch` fails with `fatal: bad object refs/Icon?`, step 2 did not run or Drive has
   already re-stamped: run it again.

4. **Open Positron** at the repo root and restore the session:
   ```r
   source(here::here("scripts", "99_session_bootstrap.R"))
   ```
   It is read-only -- no writes anywhere -- and prints a presence-and-freshness audit of the
   eleven objects the figure layer reads, the git state, and the Icon stub counts. Objects
   land in a single list `S` so nothing shadows a figure script's variables.

   For analysis work rather than figures, `options(myc.bootstrap.full = TRUE)` first; that
   also sources `00_setup_packages.R` (slow -- pulls DESeq2 and the rest).

5. **Rebuild figures only if you need to** -- see the tiers below.

---

## Getting the Claude context back

- **`claude --resume`** in the repo directory lists past sessions and restores one intact.
  The transcripts live under `~/.claude/projects/-Users-gs-G-data-MK-myc-2022-myc-mouse/`
  on the **local disk**, not on Drive, so they survive a restart untouched. The session that
  produced this runbook is `71ffe8cc-b22d-412b-9e40-f22881e23306`.
- **`claude -c`** continues the most recent session without the picker.
- `CLAUDE.md` and the memory index load automatically in any new session; the dated docs in
  `docs/` are the durable record and are all committed.
- If starting cold, paste the resume prompt at the foot of this file.

---

## Which scripts to run

**Tier 0 -- normal restart. None.** Source the bootstrap and carry on. `results/`,
`outputs/`, `paper/_freeze` and the rendered `paper/myc_mito.html` all persist on disk.

**Tier 1 -- a figure script changed.** Source that script, or all four manuscript figures at
once:
```r
source(here::here("figures", "rebuild_manuscript_figures.R"))          # writes the PDFs
options(myc.fig.nosave = TRUE)                                         # or a dry run:
source(here::here("figures", "rebuild_manuscript_figures.R"))          # draws, writes nothing
```
Re-render the document only if prose or a panel changed:
```bash
quarto render paper/myc_mito.qmd
```
Note the division of labour: **`quarto render` re-executes all eighteen figure scripts but
writes no PDFs** (its setup chunk sets `myc.fig.nosave = TRUE`, `myc_mito.qmd:28`), so it
rebuilds the *document*. `rebuild_manuscript_figures.R` is what refreshes
`outputs/figures/*.pdf`. `freeze: auto` re-executes the whole document on any source-hash
change: about 2 min changed, 29 s unchanged.

**Tier 2 -- `results/` is lost or corrupt.** Only then, and only by hand in Positron, in
this order:

| # | script | writes |
|---|--------|--------|
| 1 | `scripts/00_setup_packages.R` | -- |
| 2 | `scripts/01_load_data.R` | `count_matrix`, `coldata` |
| 3 | `scripts/03_deseq_results_qc.R` | `dds_int_run`, `interaction_results` |
| 4 | `scripts/archive_main_pipeline/02_deseq_interaction_model.R` | `combined_df_annotated{,_raw}` |
| 5 | `scripts/08_mitoPPS_analysis.R` | `mitopps_scores` |
| 6 | `scripts/32_mito_content_proxies.R` | `mito_content_proxies` |
| 7 | `scripts/40_background_vs_myc_decomposition.R` | `background_vs_myc` |
| 8 | `scripts/42_priming_arm_and_teb_substrate.R` | `priming_arm_teb` |
| 9 | `scripts/43_substrate_specificity_and_tradeoff.R` | `substrate_specificity_tradeoff` |
| 10 | `scripts/44_collapse_module_and_ownership.R` | `collapse_module_ownership` |

That closure regenerates all eleven objects the figure layer reads. Step 4 is the archived
pipeline and uses bare relative paths -- `setwd(here::here())` first. Scripts 04-07, 09-31
and 33-39 are not on the figure path; re-run one only if a specific number is needed back.

**Do not run Tier 2 to "refresh" anything.** It is hours of compute for objects that already
exist, and it overwrites values quoted verbatim in `paper/myc_mito.qmd`, in `CLAUDE.md` and
in the dated docs. The same chain is recorded in executable form at the foot of
`figures/rebuild_manuscript_figures.R`, behind an explicit opt-in.

---

## Known hazards

**The Drive `Icon\r` stub.** Google Drive for Desktop stamps a zero-byte, resource-forked
file named exactly `Icon` + CR into every folder it syncs. Anything that walks a directory
tree can trip on it. It has bitten twice:

- *Quarto*, which walks `paper/myc_mito_files/` and dies with
  `NotADirectory (os error 20)` after every chunk has already run. Handled by the
  `pre-render` / `post-render` sweeper in `paper/_quarto.yml` -- it is a **race**, not a
  state: Drive stamps the directory mid-render, so cleaning at the two ends is not enough.
- *git*, on 2026-07-27, when 272 stubs landed throughout `.git` including `refs/`. Git read
  `.git/refs/Icon\r` as a ref: `fatal: bad object refs/Icon?`. Handled by
  `./paper/clean_icon_files.sh repo`.

The durable fix is on the Drive side -- excluding this repo, or at least `.git/`, from sync.
A stub appearing in `.git/refs/` while a command is in flight is worse than a broken render.

**The stale `.RData`.** There is a 9.7 MB `.RData` in the repo root dated **2024-10-01**
(nearly two years old) and `myc_mouse.Rproj` carries `RestoreWorkspace: Default`. A front-end
that honours that would deal a year-old workspace into `globalenv()` on first launch, where
it can silently shadow correctly-loaded objects. Two ways to defuse, both the author's call:

```bash
mv .RData .RData.2024-10-01.bak     # simplest, reversible
```
or set `RestoreWorkspace: No` in `myc_mouse.Rproj` (tracked, one line).

**`combined_df_annotated{,_raw}.rds` have no writer in `scripts/`.** They come from
`scripts/archive_main_pipeline/02_deseq_interaction_model.R`, dated 2025-10-03, and are read
by `figureS1` and `figureS2`. Since `results/` is gitignored, they exist in one place only.
This is the reason the optional tar above is worth thirty seconds.

**Re-running is not free of drift.** Not every script seeds its randomness:
`21_ap6_permutation_null.R` calls `sample()` at line 125, *before* its `set.seed(1)` at line
135, and fgsea is stochastic. A re-run is not guaranteed bit-identical to the numbers now in
the manuscript. Another reason Tier 2 is a recovery path, not a refresh.

---

## Resume prompt (paste into a cold session)

> Continuing the MMTV-Myc mouse RNA-seq manuscript on branch `paper-figures`, Block B
> (publication figures + the Quarto write-up at `paper/myc_mito.qmd`).
>
> **Option A holds: you write `scripts/` + `figures/` + `paper/`, I run them in Positron;
> you verify only via scratchpad-redirected renders and never touch `outputs/` or
> `results/`.** The only untracked file,
> `docs/library_reference/gray_chea_mito_tf_shortlist.csv`, is a harmless duplicate of the
> tracked copy in `data/genesets_from_library/` -- leave it, and never `git add -A docs`.
>
> Read `CLAUDE.md` first, then `docs/2026-07-28_restart_runbook.md`. State: four manuscript
> figures built and run (`figure1`, `figure2`, `figureS1`, `figureS2`); `myc_mito.qmd` has
> no stubs left; everything committed and pushed.
>
> **Where we stopped:** I am reading the rendered Quarto and thinking about the figures, and
> will come back on that. Nothing is queued.
>
> Standing items, all author-owned (bench, not code): the Foxo3 western across both ages and
> in the PGC-1a arm; MYAZ + PGC1a RNA-seq (the only design that defines the module by
> perturbation); the METABRIC `BCL2`/`BCL2L1`/`MCL1` test; BH3 profiling of wild-type gland
> at 6W vs 12W; and the two section 0.6 citations (Emu-Myc + p53/ARF or BCL-2; beta-cell
> MycER + Bcl-xL). Open but not urgent: the paragraph 4/5 rewrite still lives only in
> section 4 of `docs/2026-07-26_introduction_alignment_and_the_question.md`, not in the
> manuscript proper.
