# c_frozen_inference.smk -- Layer C. C3 depends on B14 and has no fitting call.
rule c1_features:
    input: a10=A + "/a10/ranked_candidates.tsv.gz", fz=EXT + "/b14/frozen_M2.rds", r=EXT + "/b1_b7/resources.rds"
    output: f=A + "/c1/mva_feature_matrix.tsv.gz", manifest=A + "/c1/manifest.json"
    log: A + "/logs/c1.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/c1_features.R --config {CFG:q} --sample {wildcards.sample} --a10 {input.a10:q} --frozen {input.fz:q} --resources {input.r:q} --outdir {OUT:q}/{wildcards.sample}/c1 > {log:q} 2>&1"
rule c2_m0_m1:
    input: a10=A + "/a10/ranked_candidates.tsv.gz"
    output: m=A + "/c2/m0_m1.tsv.gz", manifest=A + "/c2/manifest.json"
    log: A + "/logs/c2.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/c2_m0_m1.R --config {CFG:q} --sample {wildcards.sample} --a10 {input.a10:q} --outdir {OUT:q}/{wildcards.sample}/c2 > {log:q} 2>&1"
rule c3_frozen_inference:
    input: f=A + "/c1/mva_feature_matrix.tsv.gz", fz=EXT + "/b14/frozen_M2.rds", fm=EXT + "/b14/freeze_manifest.json"
    output: s=A + "/c3/frozen_scores.tsv.gz", manifest=A + "/c3/manifest.json"
    log: A + "/logs/c3.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/c3_frozen_inference.R --config {CFG:q} --sample {wildcards.sample} --features {input.f:q} --frozen {input.fz:q} --freeze_manifest {input.fm:q} --outdir {OUT:q}/{wildcards.sample}/c3 > {log:q} 2>&1"
rule c4_hypotheses:
    input: a10=A + "/a10/ranked_candidates.tsv.gz", m=A + "/c2/m0_m1.tsv.gz", s=A + "/c3/frozen_scores.tsv.gz"
    output: w=A + "/c4/hypothesis_weights.tsv.gz", h=A + "/c4/final_hypothesis_set.tsv.gz", manifest=A + "/c4/manifest.json"
    log: A + "/logs/c4.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/c4_hypotheses.R --config {CFG:q} --sample {wildcards.sample} --a10 {input.a10:q} --m0m1 {input.m:q} --frozen_scores {input.s:q} --outdir {OUT:q}/{wildcards.sample}/c4 > {log:q} 2>&1"
