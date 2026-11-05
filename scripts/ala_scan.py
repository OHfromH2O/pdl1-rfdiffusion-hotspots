#!/usr/bin/env python3
"""4ZQK 복합체의 Ala 스캔 — PD-L1 접촉 잔기의 에너지 기여도 순위.

왜 돌리는가
----------
사전 등록 §1의 에피토프 핫스팟 4개(A56/A115/A121/A123)는 "PD-1과 4.5 Å 이내
접촉 + 소수성 + cleft 중심"이라는 **구조적 직관**으로 골랐습니다.
에너지 기여도로 고른 것이 아닙니다.

이 스크립트는 그 선정을 **사후 검증**합니다. 접촉 잔기를 하나씩 Ala로
치환하고 계면 에너지 변화(ddG)를 재어 순위를 냅니다.

**결과로 핫스팟을 바꾸지 않습니다.** 이미 생성된 배치의 입력을 사후에
고치면 사전 등록이 무의미해집니다. 용도는 두 가지뿐입니다.

1. 고른 4개가 에너지적으로도 상위권이었는지 확인 (→ README 한계에 반영)
2. 다음 단계 C(제3의 배치)의 입력 결정

음성 대조로 비-에피토프 핫스팟 4개(A32/A106/A109/A132)도 함께 스캔합니다.
계면에 없는 잔기이므로 ddG 가 0 근처여야 하고, 그러지 않으면 이 스캔 자체를
신뢰할 수 없습니다 — **지표의 자기검정입니다.**

주의
----
- dG_separated 는 **REU** 이며 kcal/mol 이 아닙니다. 실험 ddG 와의 상관은
  중간 수준으로 보고되어 있습니다. 순위화 보조 지표로만 씁니다.
- Ala 스캔은 사이드체인 제거의 효과를 봅니다. 원래 Ala/Gly 인 잔기는
  정의상 변화가 없으며, 이는 "기여가 없다"가 아니라 "이 방법으로 측정할 수
  없다"는 뜻입니다. **A121 은 원래 Ala 입니다** — 제외되고 그 사실을 보고합니다.
- **이완(--relax sc)을 켜면 양성 대조군이 성립합니다.** 이완 없이 4ZQK
  native 복합체를 재면 dG_separated 가 +94 REU(부호 반대)로 나옵니다.
  사이드체인만 이완하면 -32.6 REU 로 떨어지고 hbonds_int 가 5 -> 8 로
  늘어납니다. 결정 구조의 미세 충돌이 결합 상태에서 풀리지 않던 것입니다.
  다만 이완을 켜면 자기검정이 더 이상 정확히 0 이 되지 않으므로
  (전체 사이드체인이 움직입니다) 고정 임계 판정을 하지 않고 음성 대조의
  산포를 잡음 바닥으로 씁니다.
- set_pack_separated(True) 가 핵심입니다. 두 사슬을 떼어 놓은 뒤 사이드체인을
  재배치하지 않으면 분리 상태의 에너지가 실제보다 나쁘게 계산되어
  ddG 가 과대평가됩니다.

실행 (WSL2, 레포 1의 .venv-linux)
---------------------------------
    python scripts/ala_scan.py --pdb inputs/PDL1_PD1_4ZQK.pdb \
        --hotspots inputs/hotspots.json --repeats 3 --out ala_scan.csv

레포 1에서 검증된 환경: Ubuntu 22.04 (WSL2), Python 3.10, PyRosetta 2026.39.
conda base 를 비활성한 뒤 /usr/bin/python3 로 venv 를 만들어야 합니다.

**PyRosetta 라이선스는 학술·비영리 사용만 무료입니다.** 본 레포는 공개
교육 자료로서 비상업적 용도에 해당한다고 판단해 사용했습니다.
이 판단은 본 저자의 것이며 독자에게 이전되지 않습니다.
https://www.pyrosetta.org/home/licensing-pyrosetta
"""
import argparse
import csv
import json
import os
import statistics as st
import sys
import time

INTERFACE = "A_B"        # A = PD-L1 (타깃), B = PD-1 (native binder)
TARGET_CHAIN = "A"
SELFTEST_CUTOFF = 1.0    # 비계면 잔기의 |ddG| 상한 (REU)
RELAX_MODES = ("none", "sc", "all")


def init_pyrosetta():
    import pyrosetta
    # -ex1 -ex2 로 사이드체인 회전이성질체를 넓게 봅니다. 좁히면 Ala 치환 후
    # 주변 잔기가 재배치되지 못해 ddG 가 과대평가됩니다.
    pyrosetta.init("-mute all -ex1 -ex2 -use_input_sc "
                   "-ignore_unrecognized_res -ignore_zero_occupancy false")
    return pyrosetta


def make_interface_analyzer(IA):
    """생성자가 버전마다 다릅니다 (레포 1에서 겪음).
    PyRosetta 2026.x 에서 IA("A_B") 는 TypeError 를 냅니다."""
    try:
        from pyrosetta.rosetta.core.pose import DockingPartners
        dp = DockingPartners.docking_partners_from_string(INTERFACE)
        return IA(dp), "DockingPartners"
    except Exception:
        pass
    try:
        return IA(INTERFACE), "string(구버전)"
    except Exception:
        pass
    return IA(1), "jump 1"


def relax_pose(pose, sfxn, mode, repeats):
    """결합 상태를 이완시킨다.

    왜 필요한가 — 이완 없이 재면 4ZQK native 복합체의 dG_separated 가
    **+94 REU**(부호 반대)로 나옵니다. 결정 구조는 Rosetta 점수 함수에서
    보면 작은 충돌을 갖고 있고, pack_separated 로 분리 상태는 그것을 풀 수
    있지만 결합 상태는 묶여 있습니다. 그 비대칭이 양수로 나타납니다.

    mode="sc"  사이드체인만 (백본 고정). 레포 1에서 쓴 방식이고 빠릅니다.
               백본이 실험 구조이므로 건드릴 이유가 없습니다.
    mode="all" 백본·점프까지. 느리고, 계면 기하가 변하므로 H2a 좌표와
               어긋납니다. 비교용으로만.
    """
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.core.kinematics import MoveMap
    fr = FastRelax(repeats)
    fr.set_scorefxn(sfxn)
    mm = MoveMap()
    mm.set_chi(True)                  # 사이드체인은 항상 움직입니다
    mm.set_bb(mode == "all")
    mm.set_jump(mode == "all")
    fr.set_movemap(mm)
    fr.apply(pose)
    return pose


def score_interface(pose, IA_cls, sfxn):
    # 복제된 pose 는 부모의 extra score 를 물려받습니다. ia.apply 가
    # 덮어쓰지 못하면 부모(WT)의 dG_separated 를 그대로 읽어 ddG 가
    # 정확히 0 으로 나옵니다 — 실제로 겪었습니다. 먼저 비웁니다.
    try:
        from pyrosetta.rosetta.core.pose import clearPoseExtraScores
        clearPoseExtraScores(pose)
    except Exception:
        pass
    ia, _ = make_interface_analyzer(IA_cls)
    ia.set_pack_separated(True)      # 이 한 줄이 타당성을 좌우합니다
    ia.set_compute_packstat(False)   # 느리고 여기서는 쓰지 않습니다
    ia.set_scorefunction(sfxn)
    ia.apply(pose)
    d = dict(getattr(pose, "cache", None) or pose.scores)
    if "dG_separated" not in d:
        raise RuntimeError("InterfaceAnalyzerMover 가 dG_separated 를 쓰지 "
                           "않았습니다. 사슬 지정을 확인하십시오.")
    return dict(dG=d.get("dG_separated"),
                hb=d.get("hbonds_int"),
                sasa=d.get("dSASA_int"))


def pick_targets(hs, residues_arg):
    """스캔 대상과 각 잔기의 역할."""
    if residues_arg:
        tg = [int(x.strip().lstrip("Aa")) for x in residues_arg.split(",") if x.strip()]
        return tg, {r: "지정" for r in tg}
    epi_hs = {int(x[1:]) for x in hs["hotspots_epitope"]}
    ctrl_hs = {int(x[1:]) for x in hs["hotspots_control"]}
    tg = list(dict.fromkeys(list(hs["epitope_all"]) + sorted(ctrl_hs)))
    roles = {}
    for r in tg:
        if r in epi_hs:
            roles[r] = "에피토프·핫스팟"
        elif r in ctrl_hs:
            roles[r] = "비에피토프·핫스팟"
        else:
            roles[r] = "에피토프"
    return tg, roles


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", default="inputs/PDL1_PD1_4ZQK.pdb")
    p.add_argument("--hotspots", default="inputs/hotspots.json")
    p.add_argument("--repeats", type=int, default=3,
                   help="치환당 반복. 재배치가 몬테카를로이므로 중앙값을 씁니다")
    p.add_argument("--out", default="ala_scan.csv")
    p.add_argument("--residues", default="",
                   help="쉼표로 직접 지정 (기본: epitope_all + 비에피토프 핫스팟)")
    p.add_argument("--estimator", choices=("best", "median"), default="best",
                   help="best = 이완 N회 중 bound 최저를 고름 (권장). "
                        "median = 반복 dG 의 중앙값. 2026-11-04 이전의 "
                        "기본값이고 양봉에서 뒤집힙니다 — README 결과 5 의 "
                        "숫자를 재현할 때만 쓰십시오.")
    p.add_argument("--relax", choices=RELAX_MODES, default="none",
                   help="결합 상태 이완. none=기존 거동(기본), "
                        "sc=사이드체인만(권장), all=백본까지")
    p.add_argument("--relax-repeats", type=int, default=1,
                   help="FastRelax 반복. 비용이 선형으로 늘어납니다")
    a = p.parse_args()

    for f in (a.pdb, a.hotspots):
        if not os.path.isfile(f):
            sys.exit(f"{f} 가 없습니다")

    init_pyrosetta()
    from pyrosetta import pose_from_pdb, get_fa_scorefxn
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover as IA
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

    with open(a.hotspots, encoding="utf-8") as fh:
        hs = json.load(fh)
    targets, roles = pick_targets(hs, a.residues)

    sfxn = get_fa_scorefxn()
    wt = pose_from_pdb(a.pdb)
    if wt.num_chains() != 2:
        sys.exit(f"사슬이 {wt.num_chains()}개입니다. 2개여야 합니다")
    info = wt.pdb_info()

    print(f"[기준] {INTERFACE}, pack_separated=True, repeats={a.repeats}, "
          f"relax={a.relax}" + (f" x{a.relax_repeats}" if a.relax != "none" else ""))

    # parent 는 **이완하지 않은** 구조로 고정합니다. WT 기준값과 돌연변이체가
    # 똑같이 "parent 복제 -> (치환) -> 이완 1회"를 거쳐야 ddG 가 차액이 됩니다.
    #
    # 이전 판은 WT 를 먼저 이완하고 그 결과를 복제해 돌연변이시킨 뒤 또
    # 이완했습니다. 돌연변이체만 이완 2회를 받아 ddG 가 전부 음수로 끌려갔고
    # (A56 이 +3.41 -> -35.86 으로 뒤집힘), FastRelax 가 1회로 수렴하지
    # 않으므로 반복 산포가 ddG 와 같은 크기(~37 REU)로 커졌습니다.
    parent = wt
    info = parent.pdb_info()

    def prepared(mutate_at=None):
        """parent 복제 -> (치환) -> 이완 1회. WT·돌연변이체 공통 경로."""
        m = parent.clone()
        if mutate_at is not None:
            MutateResidue(mutate_at, "ALA").apply(m)
        if a.relax != "none":
            relax_pose(m, sfxn, a.relax, a.relax_repeats)
        return m

    def prepared_best(mutate_at=None):
        """이완 N회 중 **bound(결합 상태 전체 에너지)가 가장 낮은** 구조를
        고릅니다. 가장 나쁜 것도 함께 돌려줍니다 (선택이 결과를 얼마나
        바꿨는지 보이기 위해서입니다).

        왜 중앙값이 아닌가 — 이 측정에는 **서로 반대 방향의 최소화 두
        개**가 들어 있습니다. 결합 상태는 낮을수록 좋고(dG 가 더 음수로
        감), 분리 상태 재배치는 낮을수록 dG 가 **덜** 음수가 됩니다.
        그래서 dG 의 중앙값은 두 종류의 실패를 한 통에 섞습니다.

        양봉(산포 ~37 REU)의 원인은 **나쁜 이완**으로 특정되었습니다
        (README 결과 5: -70 REU 는 bound 순위 9~10위에서만 나왔고,
        IA 자체는 score-tries 5회가 모두 동일했습니다). 나쁜 이완을
        버리는 것이 올바른 선택이며, 그것이 best-of-N 입니다.
        """
        poses, bounds = [], []
        for _ in range(max(1, a.repeats)):
            m = prepared(mutate_at)
            poses.append(m)
            bounds.append(sfxn(m))      # IA 를 돌리기 전에 재야 합니다
        order = sorted(range(len(poses)), key=lambda i: bounds[i])
        return poses[order[0]], poses[order[-1]], bounds

    def measure(mutate_at=None):
        """추정량에 따라 dG 와 보조 정보를 냅니다."""
        if a.estimator == "median":
            runs = [score_interface(prepared(mutate_at), IA, sfxn)
                    for _ in range(a.repeats)]
            vals = [r["dG"] for r in runs]
            return dict(dG=st.median(vals), runs=runs, vals=vals,
                        bounds=[], dG_worst=None)
        if a.relax == "none":
            # 이완이 없으면 복제본 N 개가 전부 동일합니다. best-of-N 이
            # 무의미하므로 1 회만 재고, IA 를 헛돌리지 않습니다.
            r = score_interface(prepared(mutate_at), IA, sfxn)
            return dict(dG=r['dG'], runs=[r], vals=[r['dG']],
                        bounds=[], dG_worst=None)
        best, worst, bounds = prepared_best(mutate_at)
        rb = score_interface(best, IA, sfxn)
        rw = score_interface(worst, IA, sfxn) if a.repeats > 1 else rb
        return dict(dG=rb["dG"], runs=[rb], vals=[rb["dG"]],
                    bounds=bounds, dG_worst=rw["dG"])

    if a.relax != "none":
        # 이완이 무엇을 바꿨는지 보이지 않으면 이완을 넣은 판단을
        # 검증할 수 없습니다.
        pre = score_interface(parent.clone(), IA, sfxn)
        print(f"[이완 전] dG_separated {pre['dG']:8.2f} REU")
        t_r = time.time()
        _ = prepared()
        dt = time.time() - t_r
        print(f"[이완] {a.relax} {dt / 60:.1f}분/회 -> WT {a.repeats}회 + "
              f"치환 {len(targets)}개 x {a.repeats}회에 "
              f"약 {dt * (len(targets) + 1) * a.repeats / 60:.0f}분")

    b0 = measure()
    base_dG = b0["dG"]
    base = b0["runs"]
    if a.relax != "none":
        sign = "음수 — 양성 대조군 성립" if base_dG < 0 else "여전히 양수 — 이완으로 해결되지 않음"
        print(f"[판정] WT dG_separated {base_dG:+.2f} REU: {sign}")
    if a.estimator == "best":
        bo = ", ".join(f"{x:.0f}" for x in b0["bounds"])
        print(f"[WT] dG_separated {base_dG:8.2f} REU  (best-of-{a.repeats},"
              f" bound {bo})")
        print(f"     가장 나쁜 이완으로 재면 {b0['dG_worst']:+.2f} REU"
              f" — 차이 {abs(base_dG - b0['dG_worst']):.2f} REU")
    else:
        dg_each = ", ".join(f"{x['dG']:.2f}" for x in base)
        print(f"[WT] dG_separated 중앙값 {base_dG:8.2f} REU   (반복 {dg_each})")
    print(f"     hbonds_int {base[0]['hb']}, dSASA_int {base[0]['sasa']:.0f} A^2")
    print()

    rows = []
    t0 = time.time()
    for k, resnum in enumerate(targets, 1):
        pose_i = info.pdb2pose(TARGET_CHAIN, resnum)
        if pose_i == 0:
            print(f"  [건너뜀] A{resnum}: 구조에 없음")
            continue
        aa = parent.residue(pose_i).name3()
        if aa in ("ALA", "GLY"):
            # 잃을 사이드체인이 없습니다. 측정 불가이며 기여 없음이 아닙니다.
            rows.append(dict(res=f"A{resnum}", aa=aa, role=roles[resnum],
                             ddG="", dG_median="", n=0, spread="",
                             estimator=a.estimator, choice_gap="",
                             bound_all="", dG_all="", sasa_all="", hb_all="",
                             note=f"측정불가(원래 {aa})"))
            print(f"  [{k:2}/{len(targets)}] A{resnum} {aa}: 측정 불가 (원래 {aa})")
            continue

        m = measure(pose_i)
        runs, vals, med = m["runs"], m["vals"], m["dG"]
        ddG = med - base_dG          # 양수 = Ala 치환이 계면을 나쁘게 함 = 기여 큼
        # best 에서는 dG 반복이 1개이므로 산포를 bound 로 잽니다. 그리고
        # **이완 선택이 결과를 얼마나 바꿨는지**(최저 bound 와 최고 bound 의
        # dG 차이)를 따로 남깁니다 — 이것이 한계 11 을 직접 겨냥한 수치입니다.
        if m["bounds"]:
            spread = max(m["bounds"]) - min(m["bounds"])
            choice_gap = abs(med - m["dG_worst"])
        else:
            spread = max(vals) - min(vals)
            choice_gap = None
        # 반복별 원시값을 그대로 남깁니다. 중앙값만 보면 양봉을 놓칩니다 —
        # ~37 REU 의 상태 전환이 1/3 정도의 실행에서 일어나고 있고, 원인이
        # 아직 특정되지 않았습니다. dSASA 가 함께 뛴다면 계면 정의가
        # 바뀌는 것이고, dSASA 가 같은데 dG 만 뛴다면 분리 상태의
        # 재배치가 실패하는 것입니다.
        rows.append(dict(res=f"A{resnum}", aa=aa, role=roles[resnum],
                         ddG=round(ddG, 2), dG_median=round(med, 2),
                         n=a.repeats, spread=round(spread, 2),
                         estimator=a.estimator,
                         choice_gap=("" if choice_gap is None
                                     else round(choice_gap, 2)),
                         bound_all=";".join(f"{x:.1f}" for x in m["bounds"]),
                         dG_all=";".join(f"{v:.2f}" for v in vals),
                         sasa_all=";".join(f"{r['sasa']:.0f}" for r in runs),
                         hb_all=";".join(f"{r['hb']:.0f}" for r in runs),
                         note=""))
        # 산포가 ddG 와 같은 크기면 중앙값이 어느 봉우리에 표가 몰리는지로
        # 뒤집힙니다. 그 행은 순위로 읽을 수 없습니다.
        if choice_gap is None:
            flag = "  <- 산포가 큼, 순위 신뢰 불가" if spread > abs(ddG) else ""
            tail = f"반복 산포 {spread:.2f}"
        else:
            # choice_gap = **버려진 나쁜 분지까지의 거리**입니다. 추정값의
            # 불확실성이 아닙니다 — 처음에 그렇게 적었고 틀렸습니다.
            #
            # 실측(2026-11-04, 4ZQK 17잔기): choice_gap 이 0~2 아니면
            # 34~40 으로 **이중 분포**이고 중간이 없습니다. 그러면서도
            # 선택된 dG 는 -33.9 ~ -25.8 의 매끄러운 사다리를 이룹니다.
            # 값이 크다는 것은 나쁜 이완이 있었고 걸러냈다는 뜻이며,
            # 순위가 흔들렸다는 뜻이 아닙니다.
            #
            # 경고해야 할 경우는 반대입니다 — 두 분지가 깨끗하게 갈리지
            # 않을 때(중간값). 그러면 어느 쪽을 버려야 하는지가 모호합니다.
            flag = ("  <- 분지가 깨끗하게 갈리지 않음, 확인 필요"
                    if 5.0 < choice_gap < 20.0 else "")
            tail = f"bound 산포 {spread:.0f}, 버린 분지와의 거리 {choice_gap:.2f}"
        print(f"  [{k:2}/{len(targets)}] A{resnum} {aa}->ALA: ddG {ddG:+7.2f} REU"
              f"  (dG {med:8.2f}, {tail}){flag}")

    if not rows:
        sys.exit("계산된 잔기가 없습니다")

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # 순위 + 자기검정
    scored = sorted((r for r in rows if r["ddG"] != ""), key=lambda r: -r["ddG"])
    bar = "=" * 62
    print()
    print(bar)
    print("ddG 순위 (높을수록 계면 기여가 큼, REU)")
    print(bar)
    marks = {"에피토프·핫스팟": "  ★ 고른 것",
             "비에피토프·핫스팟": "  <- 음성 대조"}
    for i, r in enumerate(scored, 1):
        print(f"  {i:2}. {r['res']:6} {r['aa']}  {r['ddG']:+7.2f}"
              f"{marks.get(r['role'], '')}")

    picked = [r for r in scored if r["role"] == "에피토프·핫스팟"]
    ctrl = [r for r in scored if r["role"] == "비에피토프·핫스팟"]
    skipped = [r for r in rows if r["ddG"] == "" and r["role"] == "에피토프·핫스팟"]
    print()
    print("-" * 62)
    if picked:
        ranks = [scored.index(r) + 1 for r in picked]
        print(f"고른 에피토프 핫스팟의 순위: {ranks} / {len(scored)}")
    for r in skipped:
        print(f"  {r['res']} 는 순위에 없습니다 — {r['note']}")
    if ctrl:
        mx = max(abs(r["ddG"]) for r in ctrl)
        vals = ", ".join(f"{r['res']} {r['ddG']:+.2f}" for r in ctrl)
        print(f"음성 대조 (비에피토프 핫스팟): {vals}")
        print(f"  |ddG| 최대 {mx:.2f} REU  = 이 스캔의 잡음 바닥")

        if a.relax == "none":
            # 이완이 없으면 계면 밖 잔기의 내부 에너지가 dG_separated 의
            # 두 항에서 정확히 상쇄됩니다. 0 이 아니면 사슬 지정이나
            # 계면 정의가 틀렸다는 뜻이므로 고정 임계값으로 판정합니다.
            ok = "통과" if mx < SELFTEST_CUTOFF else "실패"
            print(f"  고정 임계 {SELFTEST_CUTOFF} 미만이어야 함 -> {ok}")
            if mx >= SELFTEST_CUTOFF:
                print("  비계면 잔기가 계면 에너지를 바꿉니다. 이 스캔을 신뢰할 수")
                print("  없습니다 — 사슬 지정과 pack_separated 를 확인하십시오.")
        else:
            # 이완을 켜면 set_chi(True) 가 전체 사이드체인을 움직이므로
            # 계면 밖 치환도 이완 경로를 바꾸고, 정확한 상쇄가 사라집니다.
            # 따라서 고정 임계값으로 판정할 수 없습니다. 대신 음성 대조의
            # 산포를 잡음 바닥으로 보고, 고른 핫스팟이 그것을 넘는지 봅니다.
            print(f"  relax={a.relax} 이므로 고정 임계 판정을 하지 않습니다 —")
            print("  전체 사이드체인이 움직여 계면 밖 치환도 값을 바꿉니다.")
            if picked:
                lo = min(r["ddG"] for r in picked)
                ratio = lo / mx if mx > 0 else float("inf")
                verdict = "신호 > 잡음" if lo > mx else "신호가 잡음에 묻힘"
                print(f"  고른 핫스팟 최저 {lo:+.2f} vs 잡음 {mx:.2f} "
                      f"= {ratio:.1f}배  -> {verdict}")
    print("-" * 62)
    print(f"소요 {(time.time() - t0) / 60:.1f}분  ->  {a.out}")
    print()
    print("dG_separated 는 REU 이며 kcal/mol 이 아닙니다. 순위화 보조 지표입니다.")


if __name__ == "__main__":
    main()
