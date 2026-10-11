"""R3 resolution, step 4: machine-readable spec for a corrected Fig. 1d and the corrected Sec. 4 numbers,
built only from out/r1_boot.json, out/r2_combine.json, out/r3_gaia_contrast.json and the auditor's summary.json.
Nothing in the package is modified; this is input for whoever regenerates the figure (a6_figure.py / make_figure_rev.py).
"""
import os, json
import numpy as np

A = "/Users/paulsmith/space <3/tess-jev/addenda/opus55_final_audit_tic4206066"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
r2 = json.load(open(os.path.join(OUT, "r2_combine.json")))
r3 = json.load(open(os.path.join(OUT, "r3_gaia_contrast.json")))
aud = json.load(open(A + "/03_pixels/summary.json"))["sets"]
pe = r2["per_set_errors_3.18d"]
b98 = aud["S98_clean29_B"]

spec = dict(
    note="Errors per axis = max(injection SD, jackknife SE, robust bootstrap SD); positions from audit 03 (auditor code), "
         "reproduced by the verifier's independent code within 0.3 arcsec.",
    plot=[
        dict(signal="3.18 d", sector=32, marker="filled square", product="SPOC 2-min TPF, POS_CORR-decorrelated, 8 transits",
             off=pe["S32PC"]["off_auditor"], err=pe["S32PC"]["adopted_max"]),
        dict(signal="3.18 d", sector=98, marker="filled triangle", product="TESScut 200 s, QLP flags 2048+2^29 removed, 9 transits (3995.555 excluded: no post-transit baseline)",
             off=pe["S98"]["off_auditor"], err=pe["S98"]["adopted_max"]),
        dict(signal="11.13 d", sector=98, marker="open triangle", product="TESScut 200 s, 3 clean transits",
             off=b98["real_off"], err=r2["B_S98only_auditor"]["sig"]),
    ],
    ellipses=dict(A_S32_S98=dict(mean=r2["A_eventscatter_auditor_pos"]["mean"], sig=r2["A_eventscatter_auditor_pos"]["sig"]),
                  B_S98=dict(mean=r2["B_S98only_auditor"]["mean"], sig=r2["B_S98only_auditor"]["sig"])),
    do_not_plot=dict(
        S6_both_signals="momentum dump within 0.4-3.5 h of every transit (inside the 1.4-4.7 h baseline window or in transit); position depends on recipe: "
                        "package (+8.2,-2.4), raw (-8.6,-34.5), dump-excised (+6.3,-1.3), POS_CORR (+1.0,-1.9) with jackknife SE (10.6, 32.8) arcsec",
        S32_11d="only one transit (BTJD 2197.361) passes baseline acceptance; 2175.096 lies 2.3 h after the 3.18 d transit at 2175.001",
    ),
    text_numbers=dict(
        A=dict(value="(+5 +- 3, 0 +- 3) arcsec", exact=r2["A_eventscatter_auditor_pos"], z_E=r2["A_eventscatter_auditor_pos"]["z_E"],
               p_range_over_error_models=[min(r2[k]["p"] for k in ("A_injcov_auditor", "A_eventscatter_auditor_pos", "A_bootcov_verifier_pos", "A_eventscatter_biascorr")),
                                          max(r2[k]["p"] for k in ("A_injcov_auditor", "A_eventscatter_auditor_pos", "A_bootcov_verifier_pos", "A_eventscatter_biascorr"))]),
        B=dict(value="(+1 +- 3, -1 +- 3) arcsec", exact_S98=r2["B_S98only_auditor"], exact_S32PC_S98=r2["B_S32PC+S98_auditor"]),
        capable_star_min_sigma_from_corrected=r2["capable_star_distance"]["corrected_rounded(+5+-3,0+-3)"]["min_maha"],
        gaia=dict(nearest_other_source=r3["gaia_sources_within_15as"][0],
                  depth_capable_dG_max_range=[min(x["dG_max"] for x in r3["depth_capable_limits"]), max(x["dG_max"] for x in r3["depth_capable_limits"])],
                  sep90_range_BC2019=[min(x["sep90"] for x in r3["depth_capable_limits"]), max(x["sep90"] for x in r3["depth_capable_limits"])],
                  expected_missed_chance_aligned_2p2_11as=r3["expected_chance_aligned"]["depth500_max1.0_2.2-11.0as"]["n_depth_capable_missed_by_gaia"],
                  expected_missed_chance_aligned_0_2p2as=r3["expected_chance_aligned"]["depth500_max1.0_0.0-2.2as"]["n_depth_capable_missed_by_gaia"]),
    ),
)
json.dump(spec, open(os.path.join(OUT, "r4_panel_d_spec.json"), "w"), indent=1, default=float)
print(json.dumps(spec["plot"], default=float)); print(json.dumps(spec["text_numbers"]["A"]["p_range_over_error_models"]))
print(json.dumps(spec["text_numbers"]["gaia"], default=float))
