#!/usr/bin/env python3
"""H1의 Rosetta 스테이지 — 설계된 계면의 물리 기반 점수.

사전 등록 §5의 H1은 두 지표를 씁니다.

    pAE_interaction < 10    AF2 initial guess (GPU, dl_binder_design)
    dG_separated            PyRosetta InterfaceAnalyzerMover  <- 이 스크립트

**AF2가 먼저입니다.** 여기서 쓰는 백본은 RFdiffusion 생성물이므로
dG_separated 는 "그 백본이 옳다"는 가정 위에서만 의미를 갖습니다. 그 가정을
지지하는 독립 근거(AF2 재예측)가 없으면 틀린 백본의 계면 에너지를 보고하게
됩니다. 순서를 바꾸지 마십시오 — notes/rosetta_port_plan.md 참조.

**이완이 선행합니다.** 이완 없이 재면 4ZQK native 복합체에서 +94 REU
(부호 반대)가 나옵니다. 사이드체인 이완으로 -32.7 REU 가 되어 양성 대조군이
성립했습니다. 기본값이 --relax sc 인 이유입니다.

양성 대조군 게이트
------------------
설계군을 돌리기 전에 4ZQK 전체 복합체로 지표가 작동하는지 확인합니다.

    python scripts/h1_rosetta.py --positive-control

dG_separated 가 음수가 아니면 **거기서 멈춥니다.** native 복합체가 통과하지
못하는 기준으로 설계를 판정할 수 없습니다. 이전 작업에서 AF2 의 ipTM 이
native 에서 0.12 가 나와 지표를 못 쓰게 된 것과 같은 종류의 실패이고,
같은 방식으로 미리 잡습니다.

사용
----
    # 1) 양성 대조군 먼저
    python scripts/h1_rosetta.py --positive-control

    # 2) 통과하면 설계군·대조군 (AF2 를 통과한 구조만)
    python scripts/h1_rosetta.py "af2/epi_*.pdb"  --group epitope --out h1_rosetta.csv
    python scripts/h1_rosetta.py "af2/ctrl_*.pdb" --group control --out h1_rosetta.csv

    # 3) 비교
    python scripts/h1_rosetta.py --compare --out h1_rosetta.csv

비용: 이완 0.5분/회 (사이드체인만, WSL2 2 vCPU 실측). 100건 x 3반복이면
약 2.5시간. **먼저 2건으로 재고 총량을 계산하십시오.**

주의
----
- dG_separated 는 **REU** 이며 kcal/mol 이 아닙니다. 실험 ddG 와의 상관은
  중간 수준으로 보고되어 있습니다. 군 간 비교와 순위화에만 씁니다.
- 반복 산포를 반드시 봅니다. ala_scan 에서 ~37 REU 의 양봉을 겪었고
  (README 결과 5, 한계 11) 원인이 특정되지 않았습니다. 산포가 큰 행은
  중앙값이 어느 봉우리에 표가 몰리는지로 뒤집히므로 읽을 수 없습니다.
- FastRelax 는 좌표를 바꿉니다. H2a 는 **이완 전** 좌표에서 계산했으므로
  같은 표에 놓을 때 각주를 달아야 합니다. H2a 를 재계산하지 않습니다.

**PyRosetta 라이선스는 학술·비영리 사용만 무료입니다.** 본 레포는 공개
교육 자료로서 비상업적 용도에 해당한다고 판단해 사용했습니다.
이 판단은 본 저자의 것이며 독자에게 이전되지 않습니다.
https://www.pyrosetta.org/home/licensing-pyrosetta
"""
import argparse
import csv
import glob
import os
import statistics as st
import sys
import time

# ala_scan 의 검증된 구현을 그대로 씁니다. 복사하면 고칠 곳이 둘이 됩니다.
# - make_interface_analyzer: 생성자 3단 폴백 (PyRosetta 2026.x 에서 IA('A_B') 는 TypeError)
# - score_interface:         pack_separated + extra score 비우기
# - relax_pose:              사이드체인만 이완
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ala_scan import (init_pyrosetta, relax_pose,  # noqa: E402
                      score_interface, RELAX_MODES)

POSCTRL = "inputs/PDL1_PD1_4ZQK.pdb"


def measure(pdb, IA, sfxn, relax, relax_repeats, relax_tries, score_tries):
    """한 구조의 dG_separated 를 best-of-N 으로 추정한다.

    왜 중앙값이 아닌가 — 확률 과정이 **둘**이고 방향이 반대입니다.

        dG_separated = (결합 상태 에너지) - (분리 후 재배치 에너지)

    양쪽 모두 최소화 문제이므로 각각 가장 낮은 값을 찾아야 합니다. 그런데
    분리 상태가 낮아지면 dG 는 **덜 음수**가 되고, 결합 상태가 낮아지면
    dG 는 **더 음수**가 됩니다. 둘을 섞어 중앙값을 내면 서로 상쇄되어
    어느 쪽도 최소화되지 않습니다.

    양성 대조군 10회 진단 (README 한계 11):
      - bound  는 -317 ~ -351 로 **연속** 변동 (FastRelax 의 정상 산포)
      - dG_sep 는 8회가 -32±1, 2회가 -70 인 **양봉** (간격 37 REU)
      - **-70 은 bound 가 가장 나쁜 두 실행에서 나왔습니다.** bound 순위
        9·10 위입니다. 선형 상관이 아니라 문턱이고, 이완이 충분히 안 되면
        해소되지 않은 것이 남아 dG 가 터집니다.
      - 같은 pose 에 IA 를 반복하면 **완전히 동일한 값**이 나옵니다
        (-32.74 x 5). 즉 분리 상태 재배치는 결정론적입니다.

    따라서 변동의 원천은 **이완 하나**이고, 대응은 best-of-N 입니다.

      1) relax_tries 회 이완하고 **bound 가 가장 낮은** pose 를 고른다.
         나쁜 이완(약 20%)이 걸러지고, 3회면 전부 나쁠 확률이 약 0.8% 입니다.
      2) 그 pose 를 점수화한다. score_tries 는 결정론성 확인용 카나리아이고
         (기본 2) 값이 갈리면 전제가 깨졌다는 신호입니다.

    정밀도 — 진단 10회 중 정상이던 8회의 dG 범위가 **1.7 REU** 였습니다.
    군 간 차이를 이보다 세밀하게 읽을 수 없습니다.

    경위 기록 — 처음에는 "분리 상태 재배치가 20% 확률로 실패한다"고
    진단했습니다. bound 두 값이 1.6 REU 차이인데 dG 가 37 REU 달랐기
    때문입니다. 그 두 값이 bound 순위 맨 아래 인접한 둘이라는 것을 보지
    못한 오독이었고, score_tries 가 모두 동일한 값을 내면서 드러났습니다.
    프로토콜은 결과적으로 맞았지만 이유가 틀렸습니다.

    parent 는 절대 이완하지 않습니다 — ala_scan 에서 겪은 비대칭
    (돌연변이체만 이완 2회) 을 반복하지 않기 위한 구조적 장치입니다.
    """
    from pyrosetta import pose_from_pdb
    parent = pose_from_pdb(pdb)
    if parent.num_chains() != 2:
        return None, f"사슬이 {parent.num_chains()}개 (2개여야 함)"
    # --- 1단계: 이완 best-of-N (bound 최저) ---
    best, bounds = None, []
    for _ in range(max(1, relax_tries)):
        m = parent.clone()
        if relax != "none":
            relax_pose(m, sfxn, relax, relax_repeats)
        b = sfxn(m)          # IA 를 돌리기 전에 재야 합니다 (IA 가 pose 를 건드림)
        bounds.append(b)
        if best is None or b < best[0]:
            best = (b, m)
    bound, pose = best

    # --- 2단계: 분리 상태 재배치 best-of-M (dG 가 가장 덜 음수) ---
    runs = []
    for _ in range(max(1, score_tries)):
        r = score_interface(pose.clone(), IA, sfxn)
        r["bound"] = bound
        runs.append(r)
    runs.sort(key=lambda r: r["dG"], reverse=True)   # 덜 음수가 앞
    return dict(best=runs[0], all_dG=[r["dG"] for r in runs],
                all_bound=bounds), None


def summarize(res):
    """best-of-N 결과를 한 줄로 만든다. 대표값은 중앙값이 아니라 best 입니다."""
    b, dg, bo = res["best"], res["all_dG"], res["all_bound"]
    # 분리 상태 재배치가 실패한 횟수. best 에서 5 REU 이상 더 음수인 것.
    failed = sum(1 for v in dg if v < b["dG"] - 5)
    return dict(dG=round(b["dG"], 2), bound=round(b["bound"], 2),
                hb=b["hb"], sasa=round(b["sasa"], 0),
                sep_failed=f"{failed}/{len(dg)}",
                dG_all=";".join(f"{v:.2f}" for v in dg),
                bound_all=";".join(f"{v:.2f}" for v in bo))


def compare(path):
    """CSV 에 쌓인 두 그룹을 나란히 요약한다. 계산은 하지 않는다."""
    if not os.path.isfile(path):
        sys.exit(f"{path} 가 없습니다")
    g = {}
    for r in csv.DictReader(open(path, encoding="utf-8")):
        if r["group"] in ("epitope", "control"):
            g.setdefault(r["group"], []).append(r)
    if len(g) < 2:
        sys.exit(f"두 그룹이 모두 필요합니다. 현재: {sorted(g)}")

    bar = "-" * 66
    cols = ["epitope", "control"]
    print(bar)
    print(f"H1 dG_separated 그룹 비교  ({path})")
    print(bar)
    print(f"{'':34}" + "".join(f"{c:>15}" for c in cols))

    def row(label, fn):
        print(f"  {label:32}" + "".join(f"{fn(g[c]):>15}" for c in cols))

    row("n", lambda r: len(r))
    row("dG_separated 중앙값 (REU)",
        lambda r: f"{st.median(float(x['dG']) for x in r):.2f}")
    row("  범위",
        lambda r: (f"{min(float(x['dG']) for x in r):.1f}"
                   f" ~ {max(float(x['dG']) for x in r):.1f}"))
    row("결합 상태 총 에너지 중앙값",
        lambda r: f"{st.median(float(x['bound']) for x in r):.0f}")
    row("hbonds_int 중앙값",
        lambda r: f"{st.median(float(x['hb']) for x in r):.0f}")
    row("dSASA_int 중앙값 (A^2)",
        lambda r: f"{st.median(float(x['sasa']) for x in r):.0f}")
    row("점수화가 갈린 설계 (0 이어야 정상)",
        lambda r: f"{sum(int(x['sep_failed'].split(chr(47))[0]) > 0 for x in r)}/{len(r)}")
    print()
    print("  dG_separated 는 REU 이며 kcal/mol 이 아닙니다.")
    print("  H1 의 사전 등록 기준은 pAE_interaction < 10 이고 이 값은 보조입니다.")
    print("  대표값은 중앙값이 아니라 best-of-N 입니다 (measure() docstring 참조).")
    print("  측정 정밀도는 약 +-1.7 REU 입니다. 이보다 작은 차이는 읽지 마십시오.")
    print(bar)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("pdbs", nargs="*", help="--compare / --positive-control 시 생략")
    p.add_argument("--group", default="", help="epitope | control")
    p.add_argument("--out", default="h1_rosetta.csv")
    p.add_argument("--relax", choices=RELAX_MODES, default="sc",
                   help="기본 sc. none 으로 두면 양성 대조군이 실패합니다")
    p.add_argument("--relax-repeats", type=int, default=1)
    p.add_argument("--relax-tries", type=int, default=3,
                   help="이완 횟수. bound 가 가장 낮은 것을 고릅니다")
    p.add_argument("--score-tries", type=int, default=2,
                   help="같은 pose 를 점수화하는 횟수. IA 는 결정론적이므로 "
                        "값이 갈리면 전제가 깨진 것입니다 (카나리아)")
    p.add_argument("--positive-control", action="store_true",
                   help=f"{POSCTRL} 로 지표가 작동하는지만 확인한다")
    p.add_argument("--diagnose", type=int, default=0, metavar="N",
                   help="양성 대조군을 N회 반복해 양봉의 출처를 가린다")
    p.add_argument("--compare", action="store_true")
    a = p.parse_args()

    if a.compare:
        compare(a.out)
        return

    init_pyrosetta()
    from pyrosetta import get_fa_scorefxn
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover as IA
    sfxn = get_fa_scorefxn()

    # ---------- 양봉 진단 ----------
    if a.diagnose:
        # 이완 1회 + 점수 1회를 N 쌍 돌려 둘의 분포를 **따로** 봅니다.
        # best-of-N 은 쓰지 않습니다 — 진단은 원시 분포를 보는 것입니다.
        from pyrosetta import pose_from_pdb
        print(f"[진단] {POSCTRL} x {a.diagnose}쌍, relax={a.relax}")
        print("  연속 변동과 이산 점프를 구별합니다. 산포 크기만 비교하면")
        print("  둘을 혼동합니다 — 실제로 한 번 틀린 판정을 냈습니다.")
        print()
        parent = pose_from_pdb(POSCTRL)
        rows_d = []
        for _ in range(a.diagnose):
            m = parent.clone()
            if a.relax != "none":
                relax_pose(m, sfxn, a.relax, a.relax_repeats)
            b = sfxn(m)
            r = score_interface(m, IA, sfxn)
            r["bound"] = b
            rows_d.append(r)

        print(f"  {'#':>3} {'bound (REU)':>13} {'dG_sep (REU)':>13}"
              f" {'dSASA':>8} {'hb':>4}")
        for i, r in enumerate(rows_d, 1):
            print(f"  {i:3} {r['bound']:13.2f} {r['dG']:13.2f}"
                  f" {r['sasa']:8.0f} {r['hb']:4.0f}")

        def gap(vals):
            """정렬 후 가장 큰 인접 간격과 그 위치. 양봉이면 간격이 큽니다."""
            v = sorted(vals)
            if len(v) < 3:
                return 0.0, 0
            d = [(v[i + 1] - v[i], i + 1) for i in range(len(v) - 1)]
            return max(d)

        dg = [r["dG"] for r in rows_d]
        bo = [r["bound"] for r in rows_d]
        g_dg, k_dg = gap(dg)
        g_bo, k_bo = gap(bo)
        rng_dg, rng_bo = max(dg) - min(dg), max(bo) - min(bo)
        print()
        print(f"  {'':10}{'범위':>10}{'최대 인접 간격':>16}{'간격/범위':>12}")
        print(f"  {'bound':10}{rng_bo:10.2f}{g_bo:16.2f}{g_bo / rng_bo:12.2f}")
        print(f"  {'dG_sep':10}{rng_dg:10.2f}{g_dg:16.2f}{g_dg / rng_dg:12.2f}")
        print()
        # 간격이 범위의 절반을 넘으면 두 덩어리로 갈린 것 = 양봉.
        BI = 0.5
        bi_dg, bi_bo = g_dg / rng_dg > BI, g_bo / rng_bo > BI
        if bi_dg and not bi_bo:
            print(f"  -> dG_sep 만 양봉입니다 (간격 {g_dg:.1f} REU, "
                  f"소수 {min(k_dg, len(dg) - k_dg)}/{len(dg)}건).")
            print("     bound 는 연속 변동이므로 결합 상태는 재현됩니다.")
            print("     **분리 상태 재배치의 간헐적 실패**입니다.")
            print("     -> 중앙값이 아니라 dG 가 가장 덜 음수인 것을 고르십시오")
            print("        (= 분리 상태를 가장 잘 재배치한 실행).")
        elif bi_bo:
            print("  -> bound 가 양봉입니다. FastRelax 가 서로 다른 최소점에")
            print("     떨어집니다. --relax-tries 를 늘려 최저 bound 를 고르십시오.")
        else:
            print("  -> 양쪽 모두 연속 변동입니다. 양봉이 재현되지 않았으니")
            print("     N 을 늘리거나, best-of-N 으로 충분한지 확인하십시오.")
        return

    # ---------- 양성 대조군 게이트 ----------
    if a.positive_control:
        if not os.path.isfile(POSCTRL):
            sys.exit(f"{POSCTRL} 가 없습니다")
        print(f"[양성 대조군] {POSCTRL}, relax={a.relax}, "
              f"relax-tries={a.relax_tries}, score-tries={a.score_tries}")
        t0 = time.time()
        res, err = measure(POSCTRL, IA, sfxn, a.relax, a.relax_repeats,
                           a.relax_tries, a.score_tries)
        if err:
            sys.exit(err)
        s = summarize(res)
        dt = (time.time() - t0) / 60
        print(f"  dG_separated {s['dG']:+8.2f} REU   (정밀도 약 +-1.7 REU)")
        print(f"    점수화 반복 {s['dG_all']}   <- 전부 같아야 정상")
        if s["sep_failed"] != f"0/{a.score_tries}":
            print(f"    [주의] 값이 갈렸습니다 ({s['sep_failed']}). IA 가 "
                  f"결정론적이라는 전제가 깨졌습니다.")
        print(f"  결합 상태 총 에너지 {s['bound']:+.2f} REU   <- best (최저)")
        print(f"    이완 후보 {s['bound_all']}")
        print(f"  hbonds_int {s['hb']:.0f}, dSASA_int {s['sasa']:.0f} A^2")
        print(f"  소요 {dt:.1f}분/건 -> 100건이면 약 {dt * 100 / 60:.1f}시간")
        print()
        if s["dG"] < 0:
            print("  [통과] 음수입니다. 설계군을 돌릴 수 있습니다.")
            print("  주의: PD-1 사슬 B 의 84->93 결손 8잔기(85-92)를 메우지 않은 값입니다.")
            print("  사전 등록 §4 에 따라 메운 값도 함께 보고해야 합니다.")
        else:
            print("  [실패] 양수입니다. **여기서 멈추십시오.**")
            print("  native 복합체가 통과하지 못하는 기준으로 설계를 판정할 수")
            print("  없습니다. H1 을 pAE_interaction 만으로 축소하고 그 사실을")
            print("  README 한계에 적으십시오.")
            sys.exit(1)
        return

    # ---------- 본실행 ----------
    if not a.pdbs or not a.group:
        sys.exit("pdbs 와 --group 이 필요합니다 (--compare / --positive-control 제외)")
    files = []
    for pat in a.pdbs:
        files.extend(sorted(glob.glob(pat)))
    if not files:
        sys.exit(f"일치하는 PDB 가 없습니다: {a.pdbs}")

    print(f"[{a.group}] {len(files)}건, relax={a.relax}, "
          f"relax-tries={a.relax_tries}, score-tries={a.score_tries}")
    print("  AF2 를 통과한 구조만 넣었는지 확인하십시오 — 이 스크립트는")
    print("  pAE_interaction 을 모릅니다.")
    print()

    rows, t0 = [], time.time()
    for k, f in enumerate(files, 1):
        name = os.path.splitext(os.path.basename(f))[0]
        res, err = measure(f, IA, sfxn, a.relax, a.relax_repeats,
                           a.relax_tries, a.score_tries)
        if err:
            print(f"  [{k:3}/{len(files)}] {name}: 건너뜀 — {err}")
            continue
        s = summarize(res)
        # IA 는 결정론적이어야 합니다. 갈리면 전제가 깨진 것입니다.
        nf, _ = (int(x) for x in s["sep_failed"].split("/"))
        flag = "  <- 점수화가 갈림 (전제 위반)" if nf else ""
        print(f"  [{k:3}/{len(files)}] {name}: dG {s['dG']:+8.2f} REU"
              f"  (bound {s['bound']:+.0f}, hb {s['hb']:.0f},"
              f" dSASA {s['sasa']:.0f}, 실패 {s['sep_failed']}){flag}")
        rows.append(dict(name=name, group=a.group, relax=a.relax,
                         relax_tries=a.relax_tries,
                         score_tries=a.score_tries, **s))

    if not rows:
        sys.exit("계산된 구조가 없습니다")

    new = not os.path.isfile(a.out)
    with open(a.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        if new:
            w.writeheader()
        w.writerows(rows)

    dg = [r["dG"] for r in rows]
    bad = sum(int(r["sep_failed"].split("/")[0]) > 0 for r in rows)
    print()
    print(f"  dG_separated 중앙값 {st.median(dg):+.2f} REU "
          f"(범위 {min(dg):.1f} ~ {max(dg):.1f})")
    print(f"  점수화가 갈린 설계: {bad}/{len(rows)}건 (0 이어야 정상)")
    print(f"  소요 {(time.time() - t0) / 60:.1f}분  ->  {a.out}")


if __name__ == "__main__":
    main()
