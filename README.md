Aquí tienes un **README profesional, honesto y claro** para tu repositorio **repo_cedia**, escrito exactamente como lo necesitas.
Está pensado para dejar muy buena impresión en un PI o evaluador técnico, explicando que es un workflow funcional, pragmático y 100% real-world, aunque no use Nextflow ni tecnologías modernas.

Puedes copiar–pegar directamente en `repo_cedia/README.md`.

---

# repo_cedia — Whole-Exome Sequencing (WES) Workflow

## Overview

This repository contains a **practical WES variant-calling workflow**, developed and used at *Laboratorio Biomolecular (Ecuador)* for routine whole-exome sequencing analysis.
It implements the full process **from raw FASTQ files to final variant tables ready for clinical interpretation**, typically exported to Excel or Power BI for downstream analysis.

Although the pipeline is not based on modern workflow managers (e.g., Nextflow, Snakemake), it was designed before such tools—or LLM assistance—were commonly available in the laboratory. Despite this, the workflow has proven **robust, reliable, and highly functional** in real diagnostic settings.

---

## Main Features

* **Input:** paired-end FASTQ files
* **Processing:**

  * Quality control (FastQC)
  * Alignment (BWA-MEM)
  * Sorting, indexing, duplicate marking (Samtools/Picard)
  * Coverage exploration
* **Variant calling:**

  * GATK-based workflow for SNVs and indels
  * Generation of VCF files
* **Annotation:**

  * R-based annotation modules
  * Integration of multiple databases (HPO, RefSeq, etc.)
  * Export of final interpreted variant tables
* **Output:**

  * QC reports
  * BAM/BAI
  * VCF
  * Annotated TSV/CSV for Excel
---

## Rationale

This repository reflects the **real operational constraints** of a diagnostic laboratory:

* Developed primarily in **R**, with auxiliary Bash scripts.
* Prioritized **functionality, traceability and clarity** over complex workflow automation.
* Built during a period **before LLM-based code generation existed**, where scripting pipelines manually was the standard approach.
* Structured for **daily use**, iterative updates and case-by-case analysis.

The result is a workflow that, while simple, is **practical, maintainable and validated over many real patient samples**.

---

--

## Status

This repository is intentionally **not fully refactored** or converted into a modern workflow manager, because it documents the real pipeline implemented in the laboratory.
It remains a **functional working tool** that has processed numerous WES cases in a clinical research setting.

---

## Notes

Future improvements could include:

* migration to Nextflow/Snakemake,
* Docker or Conda environments,
* modularization of R annotation scripts,
* integration with CNV-exomes workflow.

---

