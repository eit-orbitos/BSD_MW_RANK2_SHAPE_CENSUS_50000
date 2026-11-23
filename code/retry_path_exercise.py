#!/usr/bin/env python3
"""Exercise of the 384/512-bit retry path (NOT part of the study results).

The full run triggered no retry. To show that the retry ladder works, the 20 curves with the
smallest final boundary margins are re-processed with an artificially loose boundary threshold
(1e-2 instead of the predeclared 1e-20), which forces 192 -> 256 -> 384 -> 512 retries.
V0.2 labelling: the comparison with the study is between the 45-significant-digit SERIALIZED s of the
original 256-bit result and the serialized s of the rerun's final (512-bit) result. Equality of those strings
(and a zero difference of the parsed values) shows agreement to the serialized 45 digits only; it does not
establish raw 256-vs-512-bit agreement beyond that. The raw within-rerun pair diagnostics (computed on
in-memory PARI values inside the rerun) are reported separately as `within_rerun_pairs`.
Usage: python3 retry_path_exercise.py PER_CURVE_JSONL OUT_JSON"""
import sys, json
import run_study as rs
import mwshape as m

R = [json.loads(l) for l in open(sys.argv[1])]
R.sort(key=lambda r: min(float(r["final"]["margin_2b"]), float(r["final"]["margin_h"])))
sel = R[:20]
rs.BOUNDARY = "1e-2"
out = []
m.pari.set_real_precision_bits(256)   # parsing precision for serialized strings
for r in sel:
    rec = rs.process(r["source_line"] + "\n")
    s512, s256 = m.pari(rec["final"]["s"]), m.pari(r["final"]["s"])
    out.append({"label": r["label"], "margins_256": [r["final"]["margin_2b"], r["final"]["margin_h"]],
                "status_under_loose_threshold": rec["status"], "final_bits": rec["final_bits"],
                "effective_bits_final": rec["final"]["effective_bits"], "retry_events": rec["retry_events"],
                "within_rerun_pairs_raw_inmemory": [{"low": p["low_bits"], "high": p["high_bits"], "same_U": p.get("same_U"),
                                                     "s_abs": p.get("diffs", {}).get("s_abs")} for p in rec["pairs"]],
                "serialized45_s_original_256": r["final"]["s"], "serialized45_s_rerun_final": rec["final"]["s"],
                "serialized45_strings_equal": r["final"]["s"] == rec["final"]["s"],
                "abs_diff_of_parsed_serialized45_values": m.sig(abs(s512 - s256), 6)})
    print(r["label"], rec["status"], rec["final_bits"], len(rec["retry_events"]), "serialized-45 equal:", out[-1]["serialized45_strings_equal"])
json.dump({"note": "artificial threshold 1e-2 to force retries; not study results. Comparison with the study is of 45-digit serialized values only.", "cases": out},
          open(sys.argv[2], "w"), indent=1)
