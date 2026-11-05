#!/usr/bin/env python3
"""H1 의 두 지표를 태그로 조인해 한 표로 본다.

    AF2 initial guess   ->  pae_interaction, plddt_binder, binder_aligned_rmsd
    PyRosetta           ->  dG_separated, hbonds_int, dSASA_int

**사전 등록이 약속한 것과 그렇지 않은 것을 이 스크립트가 구분해 출력합니다.**

사전 등록(§4·§5)대로의 판정
---------------------------
- 주 지표는 `dG_separated` 입니다. `pae_interaction` 은 양성 대조군(native
  PD-1)이 26.06 으로 실패해 "측정 실패"로 보고하고 순위화에 쓰지 않습니다.
- 보고 형태는 **군별 중앙값**입니다. §4 가 승격 지표의 성공 정의를 적어두지
  않았으므로(= 사전 등록의 결함, README 한계 14) 이분법 판정을 하지 않습니다.
- §5 "미리 포기하는 것" 2번에 따라 **통계 검정을 하지 않습니다.** 아래의
  순위 상관(rho)은 p 값 없는 기술 통계이며 검정이 아닙니다.

사후 관찰 (사전 등록되지 않음)
------------------------------
- dSASA 층화. 임계값 500 A^2 는 **설계군 결과를 본 뒤** 분포의 빈 구간
  (249, 250, 278 | 559, 701, ...)에서 고른 것입니다. 사후 선택이므로
  보조 근거로만 읽어야 합니다. --sasa-cut 으로 바꿔 민감도를 보십시오.
- 두 지표의 순위 일치. 이것도 사후 관찰입니다.

사용
----
    python scripts/h1_join.py \
        --csv h1_rosetta.csv \
        --sc outputs/h1_af2/epi_af2.sc outputs/h1_af2/ctrl_af2.sc \
        --out h1_joined.csv

표준 라이브러리만 씁니다 — PyRosetta 환경이 아니어도 돌아갑니다.
"""
import argparse
import csv
import os
import statistics as st
import sys

GROUP_OF = {"epi": "epitope", "ctrl": "control"}
LABEL = {"epitope": "설계군 (에피토프)", "control": "대조군 (비에피토프)"}


def read_scorefile(path):
    """dl_binder_design 의 silent score 파일을 읽는다.

    헤더 이름으로 파싱합니다. 열 순서를 가정하면 조용히 어긋납니다 —
    이 파일의 첫 열은 binder_aligned_rmsd 이고 이름이 없는 열은 없습니다.
    """
    rows, hdr = [], None
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            if not ln.startswith("SCORE:"):
                continue
            f = ln.split()[1:]
            if hdr is None:
                hdr = f
                continue
            if len(f) == len(hdr):
                rows.append(dict(zip(hdr, f)))
    return rows


def num(d, k):
    try:
        return float(d[k])
    except (KeyError, TypeError, ValueError):
        return None


def spearman(xs, ys):
    """순위 상관. scipy 없이 계산합니다 (동순위는 평균 순위).

    **p 값을 내지 않습니다.** 사전 등록 §5 가 통계 검정을 포기했으므로
    기술 통계로만 씁니다.
    """
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    if len(xs) < 3:
        return None
    rx, ry = ranks(xs), ranks(ys)
    mx, my = st.mean(rx), st.mean(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    return sxy / (sxx * syy) ** 0.5


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", default="h1_rosetta.csv",
                   help="h1_rosetta.py 가 쓴 CSV")
    p.add_argument("--sc", nargs="+", required=True,
                   help="AF2 점수 파일. 파일명에 epi / ctrl 이 들어가야 합니다")
    p.add_argument("--out", default="h1_joined.csv")
    p.add_argument("--sasa-cut", type=float, default=500.0,
                   help="계면 있음/없음 경계 (A^2). 500 은 **설계군을 본 뒤** "
                        "분포의 빈 구간에서 고른 사후 임계값입니다")
    p.add_argument("--top", type=int, default=10,
                   help="순위 일치를 볼 상위 건수")
    a = p.parse_args()

    if not os.path.isfile(a.csv):
        sys.exit(f"{a.csv} 가 없습니다. 먼저 scripts/h1_rosetta.py 를 돌리십시오.")

    # --- AF2 점수: 태그 -> 행 ---
    af2 = {}
    for path in a.sc:
        if not os.path.isfile(path):
            sys.exit(f"{path} 가 없습니다")
        base = os.path.basename(path).lower()
        g = next((v for k, v in GROUP_OF.items() if k in base), None)
        if g is None:
            sys.exit(f"{path}: 파일명에서 군을 알 수 없습니다 (epi / ctrl)")
        n = 0
        for r in read_scorefile(path):
            af2[r["description"]] = dict(r, _group=g)
            n += 1
        print(f"[읽음] {path}  {n}건  -> {g}")

    # --- Rosetta CSV ---
    with open(a.csv, newline="", encoding="utf-8") as fh:
        ros = list(csv.DictReader(fh))
    print(f"[읽음] {a.csv}  {len(ros)}건")
    print()

    # --- 조인 ---
    joined, missing = [], []
    for r in ros:
        m = af2.get(r["name"])
        if m is None:
            missing.append(r["name"])
            continue
        if m["_group"] != r["group"]:
            sys.exit(f"[중단] {r['name']}: 군이 어긋납니다 "
                     f"({m['_group']} vs {r['group']}). 점수 파일을 잘못 "
                     f"짝지었습니다 — 이대로 보고하면 군이 섞입니다.")
        joined.append(dict(
            name=r["name"], group=r["group"],
            dG=float(r["dG"]), bound=float(r["bound"]),
            hb=float(r["hb"]), sasa=float(r["sasa"]),
            sep_failed=r["sep_failed"],
            pae=num(m, "pae_interaction"),
            plddt_b=num(m, "plddt_binder"),
            rmsd=num(m, "binder_aligned_rmsd"),
        ))
    if missing:
        print(f"[주의] AF2 점수에 없는 Rosetta 행 {len(missing)}건: "
              f"{missing[:5]}{' ...' if len(missing) > 5 else ''}")
        print()

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(joined[0]))
        w.writeheader()
        w.writerows(joined)

    groups = sorted({r["group"] for r in joined})
    bar = "=" * 72

    # ── 사전 등록대로의 판정 ──────────────────────────────────────────
    print(bar)
    print("사전 등록대로의 판정 — 주 지표 dG_separated, 군별 중앙값")
    print(bar)
    print(f"{'':34}" + "".join(f"{LABEL[g][:12]:>18}" for g in groups))

    def row(label, fn):
        print(f"  {label:32}" + "".join(f"{fn([r for r in joined if r['group'] == g]):>18}"
                                        for g in groups))

    row("n", lambda r: f"{len(r)}")
    row("dG_separated 중앙값 (REU)",
        lambda r: f"{st.median(x['dG'] for x in r):+.2f}")
    row("dG_separated 최저 (REU)",
        lambda r: f"{min(x['dG'] for x in r):+.2f}")
    row("hbonds_int 중앙값",
        lambda r: f"{st.median(x['hb'] for x in r):.0f}")
    row("dSASA_int 중앙값 (A^2)",
        lambda r: f"{st.median(x['sasa'] for x in r):.0f}")
    row("점수화가 갈린 설계",
        lambda r: f"{sum(int(x['sep_failed'].split('/')[0]) > 0 for x in r)}/{len(r)}")
    print()
    print("  측정 정밀도는 약 +-1.7 REU 입니다. 이보다 작은 차이는 읽지 마십시오.")
    print("  이분법 판정을 하지 않습니다 — 사전 등록 §4 가 승격 지표의 성공")
    print("  정의를 적어두지 않았고, 지금 임계값을 만들면 HARKing 입니다.")
    print()

    # ── 측정 실패로 판정된 지표 ──────────────────────────────────────
    print(bar)
    print("순위화에 쓰지 않는 지표 — pae_interaction (양성 대조군 26.06, 실패)")
    print(bar)
    row("pae_interaction 중앙값",
        lambda r: f"{st.median(x['pae'] for x in r if x['pae'] is not None):.2f}")
    row("<10 통과",
        lambda r: f"{sum(x['pae'] is not None and x['pae'] < 10 for x in r)}/{len(r)}")
    print()
    print("  관측값입니다. 판정 근거가 아닙니다 (사전 등록 §4).")
    print()

    # ── 사후 관찰 1: dSASA 층화 ──────────────────────────────────────
    print(bar)
    print(f"[사후] dSASA 층화 — 경계 {a.sasa_cut:.0f} A^2")
    print(bar)
    print("  이 경계는 **설계군 결과를 본 뒤** 분포의 빈 구간에서 골랐습니다.")
    print("  사전 등록되지 않았으므로 보조 근거로만 읽으십시오.")
    print()
    print(f"{'':34}" + "".join(f"{LABEL[g][:12]:>18}" for g in groups))
    row(f"계면 있음 (>{a.sasa_cut:.0f})",
        lambda r: f"{sum(x['sasa'] > a.sasa_cut for x in r)}/{len(r)}")
    row("  그중 dG 중앙값",
        lambda r: (f"{st.median([x['dG'] for x in r if x['sasa'] > a.sasa_cut]):+.2f}"
                   if any(x['sasa'] > a.sasa_cut for x in r) else "-"))
    row("  그중 hb 중앙값",
        lambda r: (f"{st.median([x['hb'] for x in r if x['sasa'] > a.sasa_cut]):.0f}"
                   if any(x['sasa'] > a.sasa_cut for x in r) else "-"))
    row(f"계면 없음 (<={a.sasa_cut:.0f})",
        lambda r: f"{sum(x['sasa'] <= a.sasa_cut for x in r)}/{len(r)}")
    row("  그중 dG 중앙값",
        lambda r: (f"{st.median([x['dG'] for x in r if x['sasa'] <= a.sasa_cut]):+.2f}"
                   if any(x['sasa'] <= a.sasa_cut for x in r) else "-"))
    print()
    print("  계면이 없는 설계는 계면 에너지가 0 근처입니다. '나쁘다'가 아니라")
    print("  **그 구조에 잴 계면이 없다**는 뜻이고, 전체 중앙값을 끌어당깁니다.")
    print()

    # ── 사후 관찰 2: 두 지표의 순위 일치 ─────────────────────────────
    print(bar)
    print("[사후] 두 지표가 같은 설계를 고르는가")
    print(bar)
    for g in groups:
        rs = [x for x in joined if x["group"] == g and x["pae"] is not None]
        if len(rs) < 3:
            continue
        # pae 는 낮을수록 좋고 dG 도 낮을수록 좋습니다. 둘 다 '낮을수록
        # 좋다'이므로 양의 상관이면 두 지표가 같은 방향을 가리킵니다.
        rho = spearman([x["pae"] for x in rs], [x["dG"] for x in rs])
        top_pae = {x["name"] for x in sorted(rs, key=lambda x: x["pae"])[:a.top]}
        top_dg = {x["name"] for x in sorted(rs, key=lambda x: x["dG"])[:a.top]}
        ov = len(top_pae & top_dg)
        print(f"  {LABEL[g]}")
        print(f"    Spearman rho (pae vs dG) {rho:+.2f}"
              if rho is not None else "    rho 계산 불가")
        print(f"    상위 {a.top}건 겹침 {ov}/{a.top}"
              f"  (무작위 기대 {a.top * a.top / len(rs):.1f})")
        print()
    print("  p 값을 내지 않습니다 — 사전 등록 §5 가 통계 검정을 포기했습니다.")
    print("  rho 는 기술 통계이고, 겹침 수도 기대값과 나란히만 봅니다.")
    print()

    # ── 상위 목록 ────────────────────────────────────────────────────
    print(bar)
    print(f"dG_separated 상위 {a.top}건 (군 표시)")
    print(bar)
    print(f"  {'설계':34}{'dG':>9}{'hb':>5}{'dSASA':>8}{'pae':>8}"
          f"{'pLDDT_b':>9}{'rmsd':>7}")
    for x in sorted(joined, key=lambda x: x["dG"])[:a.top]:
        tag = x["name"].replace("_dldesign_0_cycle1_af2pred", "")
        mark = " *" if x["group"] == "epitope" else "  "
        pae = f"{x['pae']:.2f}" if x["pae"] is not None else "-"
        pb = f"{x['plddt_b']:.1f}" if x["plddt_b"] is not None else "-"
        rm = f"{x['rmsd']:.2f}" if x["rmsd"] is not None else "-"
        print(f"{mark}{tag:34}{x['dG']:+9.2f}{x['hb']:5.0f}{x['sasa']:8.0f}"
              f"{pae:>8}{pb:>9}{rm:>7}")
    print()
    print("  * = 설계군(에피토프)")
    print(bar)
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
