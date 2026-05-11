```
▄ ▀     ╔════╗      ╔═════╗         ╔════╗      ╔════╗ 
    ██████╗* ║   ▓▓▓▓▓╗ * ╚═▒▒▒╗   ▒▒▒╗  ║ ░░░░░░░╗ *║    
 ▌ ██╔════╝  ╚╗ ▓▓╔══▓▓╗    ▒▒▒▒╗ ▒▒▒▒║  ║ ░░╔════╝ *║    
   ██║  ███╗  ║ ▓▓▓▓▓▓▓║    ▒▒╔▒▒▒▒╔▒▒║  ║ ░░░░░╗  ╔═╝    
╔══██║   ██║  ║ ▓▓╔══▓▓║  * ▒▒║╚▒▒╔╝▒▒║  ╚═░░╔══╝  ╚═╗   
║* ╚██████╔╝  ╚═▓▓║  ▓▓║    ▒▒║ ╚═╝ ▒▒║  * ░░░░░░░╗  ║    
║   ╚═════╝  ** ╚═╝  ╚═╝    ╚═╝     ╚═╝    ╚══════╝  ║ 
╚════════════════════════════════════════════════════╝
````
**Genomes Analysis Made Easy**

This repository contains a scalable, user-friendly pipeline for standardised, parallel analysis of genomes being developed as part of [BGE+](https://biodiversitygenomics.eu/2026/02/16/bgeplus-selected-for-funding/) and building upon lessons learned from [ERGA](https://www.erga-biodiversity.eu/team-1/sac---sequencing-and-assembly-committee) and [GEP2](https://github.com/diegomics/GEP2).

Data is entered via a simple table, and configuration is managed through a tidy control panel. GAME uses a modern [Snakemake](https://snakemake.readthedocs.io) version with [containers](https://apptainer.org) and can run on a server/cluster (SLURM) or a local computer.

**Please cite:** *in preparation*

---

## Requirements

- `Conda`*
- `Apptainer`

*(or you could have the packages listed in the install.yml, in addition to Apptainer/Singularity, installed in your PATH)

---

## GAME can:
```
• download assemblies & reads (or use the ones in your local storage)
• trim/filter/qc reads (paired-end, 10x, HiFi, ONT), plus QC reports
• produce a quick QC report to check the reference assembly
• de novo masking with RepeatModeler/Masker (or use your bed masking or extract softmasked from reference)
• map reads (paired-end, 10x, HiFi, ONT) with dynamic merging, plus QC reports
• variant calling with DeepVariant or GATK, plus Joint Genotyping with GLnexus
• apply tags for flexible downstream soft filtering, plus QC reports
• check relatedness and LD prunning tagging, plus QC reports
• analyse heterozygosity and ROH (currently adding functionalities. Coming with next version bump!)
• generate PSMC analysis, plots and diagnostics (currently adding other demgraphic inference methods, like GONE, etc. Coming with further version bump!)
• more analysis coming with further version bumps! (Fst, PCA, admixture, structure...)
```

<Br />

### Please visit the [GAME Wiki](https://github.com/diegomics/GAME/wiki) for detailed information on how to get, set up and run this cool pipeline!

<Br />
