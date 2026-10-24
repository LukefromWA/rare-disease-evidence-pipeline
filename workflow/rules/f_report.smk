# f_report.smk -- LaTeX reporting. Assets are machine-generated from manifests/outputs; a LaTeX failure fails only the report target.
rule report_assets:
    input: a10=OUT + "/{sample}/a10/ranked_candidates.tsv.gz", c4=OUT + "/{sample}/c4/final_hypothesis_set.tsv.gz", d7=OUT + "/{sample}/d7/final_report.tsv"
    output: done="reports/generated/.assets_{sample}"          # values.tex and generated/*.tex are written as a side effect of the same node
    log: "reports/output/{sample}_assets.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/make_report_assets.R --config {CFG:q} --sample {wildcards.sample} --outroot {OUT:q}/{wildcards.sample} --reports reports > {log:q} 2>&1 && touch {output.done:q}"
rule report_track1:
    input: tex="reports/track1_report.tex", done="reports/generated/.assets_{sample}", manifest=OUT + "/{sample}/a10/manifest.json"
    output: pdf="reports/output/{sample}_track1_report.pdf"
    log: "reports/output/{sample}_track1_report.log"
    shell: "bash reports/build_report.sh {input.tex:q} {output.pdf:q} > {log:q} 2>&1"
