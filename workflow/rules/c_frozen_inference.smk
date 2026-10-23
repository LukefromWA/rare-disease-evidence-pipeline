# c_frozen_inference.smk -- Layer C. C3 depends on B14 and has no fitting call.
rule c1_features:
    input: a10=A + "/a10/ranked_candidates.tsv.gz", fz=EXT + "/b14/frozen_M2.rds", r=EXT + "/b1_b7/resources.rds"
    output: f=A + "/c1/mva_feature_matrix.tsv.gz", manifest=A + "/c1/manifest.json"
    log: A + "/logs/c1.log"
    shell: "{RSCRIPT} {SCRIPTS}/c1_features.R --config {CFG} --sample {wildcards.sample} --a10 {input.a10} --frozen {input.fz} --resources {input.r} --outdir {OUT}/{wildcards.sample}/c1 > {log} 2>&1"
rule c2_m0_m1:
    input: a10=A + "/a10/ranked_candidates.tsv.gz"
    output: m=A + "/c2/m0_m1.tsv.gz", manifest=A + "/c2/manifest.json"
    log: A + "/logs/c2.log"
    shell: "{RSCRIPT} {SCRIPTS}/c2_m0_m1.R --config {CFG} --sample {wildcards.sample} --a10 {input.a10} --outdir {OUT}/{wildcards.sample}/c2 > {log} 2>&1"
rule c3_frozen_inference:
    input: f=A + "/c1/mva_feature_matrix.tsv.gz", fz=EXT + "/b14/frozen_M2.rds", fm=EXT + "/b14/freeze_manifest.json"
    output: s=A + "/c3/frozen_scores.tsv.gz", manifest=A + "/c3/manifest.json"
    log: A + "/logs/c3.log"
    shell: "{RSCRIPT} {SCRIPTS}/c3_frozen_inference.R --config {CFG} --sample {wildcards.sample} --features {input.f} --frozen {input.fz} --freeze_manifest {input.fm} --outdir {OUT}/{wildcards.sample}/c3 > {log} 2>&1"
rule c4_hypotheses:
    input: a10=A + "/a10/ranked_candidates.tsv.gz", m=A + "/c2/m0_m1.tsv.gz", s=A + "/c3/frozen_scores.tsv.gz"
    output: w=A + "/c4/hypothesis_weights.tsv.gz", h=A + "/c4/final_hypothesis_set.tsv.gz", manifest=A + "/c4/manifest.json"
    log: A + "/logs/c4.log"
    shell: "{RSCRIPT} {SCRIPTS}/c4_hypotheses.R --config {CFG} --sample {wildcards.sample} --a10 {input.a10} --m0m1 {input.m} --frozen_scores {input.s} --outdir {OUT}/{wildcards.sample}/c4 > {log} 2>&1"
