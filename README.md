# Ampliseq @ Rieux — multi-region CLI wrapper

A bash wrapper that runs [nf-core/ampliseq](https://nf-co.re/ampliseq) across
every region of a **multi-region 16S/ITS amplicon panel** on the Rieux HPC
cluster, optionally runs the multi-region reconstruction with
[q2-sidle](https://q2-sidle.readthedocs.io), and consolidates everything into a
single `final_reports/` directory.

Built for the **QIAseq 16S/ITS Pro Screening Panel (96)** — six 16S regions
(V1V2, V2V3, V3V4, V4V5, V5V7, V7V9) plus ITS1.

What is panel-specific and what is not:

- The **primers** in `assets/primers_painel.tsv` are that kit's, recovered from
  sequencing data with `bin/discover_primers.py` (QIAGEN does not publish them).
  Another panel supplies its own table — or runs that script on its own reads.
- The **truncation lengths** in `assets/parametros_regioes.tsv` were derived
  from one specific run. They depend on read length and quality, so they are a
  starting point, not a setting: re-derive them for your run with
  `bin/perfil_qualidade.py`. Getting this wrong is expensive — an earlier
  version of that table capped V1V2 below the real community median and
  silently destroyed 98% of its merges.
- The **control name** defaults to the kit's Smart Control (`^[Ss]mart`) and is
  a flag everywhere (`--controles`), because what counts as a control belongs to
  the experiment, not to the tool.
- The **region names** come from the primers table, not from a naming
  convention. A panel whose regions are called `region1..region5` — the example
  in ampliseq's own docs — works unchanged.

`final_reports/` is the contract with the companion R package
[`aspp`](https://github.com/thiagoparentefiocruz/aspp), which reads those TSVs
and produces the figures and the downstream analysis. Same split as
`rieux_bacass` → `bpp`: the wrapper extracts metrics on the cluster with the
standard library only, the R package does everything that needs R.

## Installation

```bash
git clone https://github.com/thiagoparentefiocruz/rieux_ampliseq.git ~/rieux_ampliseq
chmod +x ~/rieux_ampliseq/rieux_ampliseq.sh ~/rieux_ampliseq/bin/*
echo 'export PATH="$HOME/rieux_ampliseq:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

Then tell it, once, where your installation lives — containers, `NXF_HOME` and
the reference databases:

```bash
echo 'RIEUX_PIPELINE_BASE=/path/to/your/pipeline' > ~/.rieux_ampliseq.conf
```

No path in this repository is tied to any particular account. The wrapper
sources `bin/ambiente.sh` on its own if the reference-database variables are not
already in the environment, and `ambiente.sh` refuses to run rather than guess
that directory.

## Invocation

The Nextflow driver runs on the **login node**, inside a `screen`. It sits idle
waiting on SLURM and submits the tasks; it needs no allocation of its own.

```bash
screen -S fabio
rieux_ampliseq.sh --name fabio --input split/samplesheets/fabio
# Ctrl-A then D
```

Regions run **in sequence**, on purpose: two regions of the same dataset would
only compete for the same queue. Two **datasets** in parallel is a different
matter — different partitions, and there the gain is real:

```bash
rieux_ampliseq.sh --name fabio    --input .../fabio                        # cpu
rieux_ampliseq.sh --name patricia --input .../patricia \
                  --partition fat --work-dir exec/patricia
```

`--work-dir` is not optional in that case: Nextflow keeps `.nextflow.log`,
`.nextflow/history` and the `-resume` cache in the launch directory, and two
drivers in the same directory shred each other's history.

Every run is idempotent — each region gets `-resume` and its own `outdir`, so
re-running after a failure picks up where it stopped.

## Key flags

| flag | meaning |
|---|---|
| `--name NAME` | label for the dataset; goes into output paths and logs |
| `--input DIR` | directory holding one `samplesheet_<REGION>.tsv` per region |
| `--outdir DIR` | output root (default `./resultados/<NAME>`) |
| `--partition NAME` | SLURM partition (default `cpu`) |
| `--regions "A B"` | only these regions |
| `--work-dir DIR` | Nextflow launch directory |
| `--multiregion FILE` | enable the Sidle branch (see below) |
| `--sidle-input FILE` | samplesheet of the **undivided** reads |
| `--dry-run` | print the commands, run nothing |

`--help` lists the rest.

## The Sidle branch

`--multiregion` takes the `regions_multiregion.tsv` produced by
`bin/validar_regioes_sidle.py`, which is a step you should not skip. That file
carries, per region, the primer pair **and** a `region_length`, and ampliseq
trims every sequence to that length and **discards anything shorter**.

The validator does in-silico PCR against the reference database and compares
the resulting length distribution with the observed ASV lengths, so both the
primer boundaries and the cut are chosen from data:

```bash
bin/validar_regioes_sidle.py \
    --primers assets/primers_painel.tsv \
    --ref "$DB_SILVA_GENERO" \
    --asv resultados/fabio/final_reports/comprimento_asv.tsv \
    --out regions_multiregion.tsv
```

The Sidle branch is fed the **undivided** reads, not the per-region split: it
does its own routing, running cutadapt once per region with those primers.
Feeding it already-split reads would trim the primers twice and shift the
boundaries — exactly what the validator exists to prevent.

Sidle is reference-bound: an ASV that matches nothing in the reference simply
disappears. That is why the per-region branch is not an alternative to it but
its denominator — without it you cannot know what the reconstruction dropped.

## Output structure

```
<outdir>/<REGION>/          one ampliseq run per region
<outdir>/sidle/             the multi-region branch, if enabled
<outdir>/logs/
<outdir>/final_reports/     the consolidated tables
```

`final_reports/` holds six TSVs:

| file | one row per |
|---|---|
| `reads_por_regiao.tsv` | sample × region, through every DADA2 step |
| `classificacao.tsv` | region × rank — how far classification got |
| `abundancia.tsv` | taxon × region, aggregated, controls separated |
| `comprimento_asv.tsv` | ASV length histogram per region |
| `abundancia_por_amostra.tsv` | taxon × sample |
| `prevalencia.tsv` | taxon × region — in how many samples it occurs |

The last two answer a question the aggregated tables cannot: a taxon at 42% may
be in every sample or piled into one, and those call for opposite decisions.

Relative abundance is computed **within region** — comparing raw counts across
regions would compare primer efficiency, not composition. Chloroplast is
filtered by **Order** (it sits inside Cyanobacteriota in SILVA, so filtering by
phylum would erase real cyanobacteria) and Mitochondria by **Family** (it sits
in Rickettsiales, the same order as Anaplasmataceae — filtering by order would
erase the target of anyone looking for *Anaplasma*).

## Partition handling

There is one config file, `conf/rieux.config`. The partition, the queue size
and the `process_high` memory ceiling come from environment variables that the
wrapper sets from `--partition`; called directly from `nextflow`, the file falls
back to the `cpu` defaults. `DADA2_TAXONOMY`/`DADA2_ADDSPECIES` always go to the
fat partition — with SILVA they are the only genuinely memory-hungry step.

## Utilities in `bin/`

| script | what it does |
|---|---|
| `ambiente.sh` | environment for the Nextflow driver |
| `split_regioes.sh` | routes raw reads into per-region samplesheets |
| `discover_primers.py` | recovers the panel's primers from the FASTQs |
| `validar_regioes_sidle.py` | in-silico PCR; writes `regions_multiregion.tsv` |
| `fazer_samplesheet.py` | builds the undivided-reads samplesheet for the Sidle branch |
| `preparar_unite.py` | formats UNITE into the two FASTA files DADA2 needs |
| `coletar_metricas.py` | builds `final_reports/` |
| `avaliar_execucao.py` | per-run QC report |
| `perfil_qualidade.py` | quality profiles, for choosing truncation |

### One gotcha worth stating out loud

`--dada_ref_tax_custom` **does not** run the `fmtscript` that
`--dada_ref_taxonomy` runs. Passing a database by hand means passing it already
formatted — and nothing warns you if you pass it raw. `preparar_unite.py` exists
because of this: it reproduces ampliseq's own `taxref_reformat_unite.sh` and
emits the **two** files (`assignTaxonomy` and `addSpecies`) that the official
path uses.

## License

MIT.
