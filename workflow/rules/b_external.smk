# b_external.smk -- Layer B (external only). Order: B1-B7 -> B8 -> B10 -> B9 -> B11 -> B12 -> B12.5 -> B13 -> B14. No private input anywhere.
rule b1_b7_prepare:
    output: c=EXT + "/b1_b7/cases.tsv.gz", r=EXT + "/b1_b7/resources.rds", manifest=EXT + "/b1_b7/manifest.json"
    log: EXT + "/logs/b1_b7.log"
    shell: "{RSCRIPT} {SCRIPTS}/b1_b7_prepare_external.R --config {CFG} --outdir {EXT}/b1_b7 > {log} 2>&1"
rule b8_candidates:
    input: c=EXT + "/b1_b7/cases.tsv.gz", r=EXT + "/b1_b7/resources.rds"
    output: x=EXT + "/b8/candidate_sets.tsv.gz", manifest=EXT + "/b8/manifest.json"
    log: EXT + "/logs/b8.log"
    shell: "{RSCRIPT} {SCRIPTS}/b8_candidate_sets.R --config {CFG} --cases {input.c} --resources {input.r} --outdir {EXT}/b8 > {log} 2>&1"
rule b10_splits:
    input: x=EXT + "/b8/candidate_sets.tsv.gz"
    output: s=EXT + "/b10/split_manifest.tsv", manifest=EXT + "/b10/manifest.json"
    log: EXT + "/logs/b10.log"
    shell: "{RSCRIPT} {SCRIPTS}/b10_splits.R --config {CFG} --candidates {input.x} --outdir {EXT}/b10 > {log} 2>&1"
rule b9_features:
    input: c=EXT + "/b1_b7/cases.tsv.gz", x=EXT + "/b8/candidate_sets.tsv.gz", s=EXT + "/b10/split_manifest.tsv", r=EXT + "/b1_b7/resources.rds"
    output: f=EXT + "/b9/features.tsv.gz", st=EXT + "/b9/feature_state.rds", manifest=EXT + "/b9/manifest.json"
    log: EXT + "/logs/b9.log"
    shell: "{RSCRIPT} {SCRIPTS}/b9_features.R --config {CFG} --cases {input.c} --candidates {input.x} --splits {input.s} --resources {input.r} --outdir {EXT}/b9 > {log} 2>&1"
rule b11_fit_m2:
    input: f=EXT + "/b9/features.tsv.gz"
    output: m=EXT + "/b11/m2_fits.rds", c=EXT + "/b11/m2_coefficients.tsv", manifest=EXT + "/b11/manifest.json"
    log: EXT + "/logs/b11.log"
    shell: "{RSCRIPT} {SCRIPTS}/b11_fit_m2.R --config {CFG} --features {input.f} --outdir {EXT}/b11 > {log} 2>&1"
rule b12_compare:
    input: f=EXT + "/b9/features.tsv.gz", m=EXT + "/b11/m2_fits.rds"
    output: r=EXT + "/b12/dependence_report.tsv", s=EXT + "/b12/selected_M2.rds", manifest=EXT + "/b12/manifest.json"
    log: EXT + "/logs/b12.log"
    shell: "{RSCRIPT} {SCRIPTS}/b12_compare_dependence.R --config {CFG} --features {input.f} --fits {input.m} --outdir {EXT}/b12 > {log} 2>&1"
rule b12_5_conformal:
    input: f=EXT + "/b9/features.tsv.gz", s=EXT + "/b12/selected_M2.rds"
    output: c=EXT + "/b12_5/conformal.json", manifest=EXT + "/b12_5/manifest.json"
    log: EXT + "/logs/b12_5.log"
    shell: "{RSCRIPT} {SCRIPTS}/b12_5_conformal.R --config {CFG} --features {input.f} --selected {input.s} --outdir {EXT}/b12_5 > {log} 2>&1"
rule b13_hierarchy:
    input: s=EXT + "/b12/selected_M2.rds"
    output: h=EXT + "/b13/hierarchy_report.tsv", manifest=EXT + "/b13/manifest.json"
    log: EXT + "/logs/b13.log"
    shell: "{RSCRIPT} {SCRIPTS}/b13_hierarchy.R --config {CFG} --outdir {EXT}/b13 > {log} 2>&1"
rule b14_freeze:
    input: s=EXT + "/b12/selected_M2.rds", st=EXT + "/b9/feature_state.rds", c=EXT + "/b12_5/conformal.json", h=EXT + "/b13/hierarchy_report.tsv"
    output: m=EXT + "/b14/frozen_M2.rds", fm=EXT + "/b14/freeze_manifest.json", sums=EXT + "/b14/checksums.sha256", manifest=EXT + "/b14/manifest.json"
    log: EXT + "/logs/b14.log"
    shell: "{RSCRIPT} {SCRIPTS}/b14_freeze.R --config {CFG} --selected {input.s} --state {input.st} --conformal {input.c} --outdir {EXT}/b14 > {log} 2>&1"
