#!/usr/bin/env python3
"""H2 — 설계된 binder가 지정한 핫스팟을 실제로 접촉하는가.

사전 등록 §5의 H2를 계산합니다. H1(pAE_interaction)보다 직접적인
"지정 효과"의 측정이며, RFdiffusion 출력 PDB만으로 계산됩니다 —
ProteinMPNN·AF2 없이 바로 돌 수 있습니다.

정의 — 두 가지 모드
------------------
**H2a (백본 프록시, 기본).** CA–CA 8.0 Å 이내.
RFdiffusion 출력에는 **사이드체인이 없습니다**(N/CA/C/O만). 계면을 건너
닿는 것은 사이드체인이므로, 백본만으로 4.5 Å을 재면 거의 0이 나옵니다.
CA–CA 8 Å은 백본만 있을 때 쓰는 관용적 접촉 정의입니다.

**H2b (원 기준).** 중원자 4.5 Å 이내. 사이드체인이 생긴 뒤
(ProteinMPNN·FastRelax·AF2 복합체)에만 쓸 수 있습니다.

사이드체인이 있는지 자동으로 판별해 모드를 고릅니다.
`--atoms`/`--cutoff`로 강제할 수 있습니다.

- 성공: 지정 핫스팟 4개 중 2개 이상 접촉 (두 모드 공통)

잔기 번호 주의
--------------
RFdiffusion 출력의 타깃 사슬은 입력 번호(18–132)를 그대로 쓰지 않을 수
있습니다. `.trb`에 원본↔출력 잔기 대응이 들어 있으므로 그것을 우선 읽고,
없으면 순서 기반으로 매핑한 뒤 그 사실을 보고합니다.

사용
----
    python scripts/h2_hotspot_contacts.py outputs/pdl1_epi_*.pdb \
        --hotspots A56,A115,A121,A123 --group epitope
    python scripts/h2_hotspot_contacts.py outputs/pdl1_ctrl_*.pdb \
        --hotspots A32,A106,A109,A132 --group control
두 결과를 같은 CSV에 append 하려면 `--out h2_contacts.csv`를 같게 주십시오.
"""
import argparse, csv, glob, os, pickle, sys
from collections import defaultdict

MIN_HIT = 2           # 4개 중 2개 이상이면 성공. 사전 등록 §5에 고정

# 사이드체인 유무에 따른 기준. 사전 등록 §5 참조.
MODES = {
    "heavy4.5": dict(cutoff=4.5, ca_only=False, label="H2b 중원자 4.5 A"),
    "ca8":      dict(cutoff=8.0, ca_only=True,  label="H2a CA-CA 8.0 A"),
}
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def read_pdb(path, ca_only):
    """사슬별 {잔기번호: [(x,y,z), ...]}. (chains, 사이드체인_있음)"""
    chains = defaultdict(lambda: defaultdict(list))
    has_sc = False
    with open(path) as fh:
        for l in fh:
            if not l.startswith("ATOM"):
                continue
            an = l[12:16].strip()
            if l[76:78].strip() == "H" or an.startswith("H"):
                continue          # 중원자만
            if an not in BACKBONE:
                has_sc = True
            if ca_only and an != "CA":
                continue
            chains[l[21]][int(l[22:26])].append(
                (float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return chains, has_sc


def target_numbering(pdb_path, target_chain, target_resnums, first=18):
    """출력 번호 → 원본 번호 매핑. (dict, 출처 문자열)"""
    trb = os.path.splitext(pdb_path)[0] + ".trb"
    if os.path.isfile(trb):
        try:
            with open(trb, "rb") as fh:
                d = pickle.load(fh)
            ref = d.get("complex_con_ref_pdb_idx")
            hal = d.get("complex_con_hal_pdb_idx")
            if ref and hal and len(ref) == len(hal):
                m = {h[1]: r[1] for r, h in zip(ref, hal)
                     if h[0] == target_chain}
                if m:
                    return m, "trb"
        except Exception as e:
            print(f"  [경고] {os.path.basename(trb)} 읽기 실패: {e}",
                  file=sys.stderr)
    # 폴백 — .trb 가 없을 때.
    #
    # 2026-11-05 수정. 옛 조건은 lo >= 18 이면 "이미 원본 번호"로 보고
    # 그대로 썼습니다. AF2 initial guess 출력은 타깃 사슬이 binder 뒤로
    # 밀려 79-193 처럼 나오므로(binder 78잔기 + 타깃 115잔기) 그 조건을
    # 통과하고, **79번을 PD-L1 79번으로 착각합니다.** 실제로는 18번입니다.
    # 조용히 틀린 잔기를 세게 되므로 조건을 lo == 18 로 좁혔습니다.
    lo = min(target_resnums)
    if lo == first:
        return {n: n for n in target_resnums}, "as-is"
    m = {n: first + i for i, n in enumerate(sorted(target_resnums))}
    print(f"  [번호] 타깃 사슬이 {lo}-{max(target_resnums)} 입니다. "
          f"순서 기반으로 {first}- 에 맞춥니다 (오프셋 {lo - first}). "
          f"타깃이 연속이라는 가정입니다 — 결손이 있으면 틀립니다.",
          file=sys.stderr)
    return m, f"positional({first}-)"


def contacts(binder_atoms, target_res, cutoff):
    """cutoff 이내인 타깃 잔기 번호 집합. 네이티브 파이썬만 사용."""
    c2 = cutoff * cutoff
    hits = set()
    for rid, atoms in target_res.items():
        for bx, by, bz in binder_atoms:
            for ax, ay, az in atoms:
                if (ax-bx)**2 + (ay-by)**2 + (az-bz)**2 <= c2:
                    hits.add(rid)
                    break
            if rid in hits:
                break
    return hits


def compare(path):
    """CSV에 쌓인 두 그룹을 나란히 요약한다. 계산은 하지 않는다."""
    import statistics as st
    if not os.path.isfile(path):
        sys.exit(f"{path} 가 없습니다. 먼저 각 그룹을 계산하십시오.")
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    groups = []
    for g in dict.fromkeys(r["group"] for r in rows):   # 입력 순서 유지
        groups.append((g, [r for r in rows if r["group"] == g]))
    if len(groups) < 2:
        sys.exit(f"그룹이 {len(groups)}개뿐입니다. 두 그룹을 모두 계산하십시오.")

    def col(r, k):
        return [int(x[k]) for x in r]

    bar = "-" * 58
    print(f"\n{bar}\nH2a 그룹 비교  ({path})\n{bar}")
    print(f"{'':28}" + "".join(f"{g:>15}" for g, _ in groups))

    def line(label, fn):
        print(f"{label:28}" + "".join(f"{fn(r):>15}" for _, r in groups))

    line("n", lambda r: len(r))
    line("핫스팟 2개 이상 (= H2a)",
         lambda r: f"{sum(1 for x in col(r,'success') if x)}/{len(r)}"
                   f" {sum(col(r,'success'))/len(r):.0%}")
    line("  3개 이상",
         lambda r: f"{sum(1 for h in col(r,'n_hotspot_hit') if h>=3)}/{len(r)}")
    line("핫스팟 접촉 수 중앙값",
         lambda r: f"{st.median(col(r,'n_hotspot_hit')):.0f}")
    line("계면 잔기 수 중앙값",
         lambda r: f"{st.median(col(r,'n_contact_total')):.0f}")
    line("  범위",
         lambda r: f"{min(col(r,'n_contact_total'))}-{max(col(r,'n_contact_total'))}")
    line("binder 길이 중앙값",
         lambda r: f"{st.median(col(r,'binder_len')):.0f}")

    print(f"\n접촉 수 분포 (0 / 1 / 2 / 3 / 4)")
    for g, r in groups:
        h = col(r, "n_hotspot_hit")
        print(f"  {g:12} " + " / ".join(str(h.count(k)) for k in range(5)))

    print(f"\n핫스팟별 접촉률")
    for g, r in groups:
        labels = r[0].get("hotspots", "").split(";")
        cols = [f"hit{i+1}" for i in range(len(labels))]
        print(f"  {g:12} " + "  ".join(
            f"{lab} {sum(int(x[c]) for x in r)}/{len(r)}"
            for lab, c in zip(labels, cols)))

    # 사전 등록 §5 판정
    rates = [sum(col(r, "success")) / len(r) for _, r in groups]
    print(f"\n{bar}")
    if len(rates) == 2 and rates[1] > 0:
        ratio = rates[0] / rates[1]
        ok = "통과" if ratio >= 2.0 else "불통과"
        print(f"사전 등록 §5 H2: {groups[0][0]} / {groups[1][0]} = "
              f"{ratio:.2f}배  (기준 2.0배 이상)  -> {ok}")
    print("계면 잔기 수 차이는 사전 등록된 기준이 아닙니다 (사후 관찰).")
    print(f"{bar}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("pdbs", nargs="*", help="--compare 시에는 생략")
    p.add_argument("--hotspots", default="",
                   help="예: A56,A115,A121,A123")
    p.add_argument("--group", default="", help="epitope | control")
    p.add_argument("--out", default="h2_contacts.csv")
    p.add_argument("--target-chain", default=None,
                   help="기본: 잔기 수가 많은 사슬을 타깃으로 본다")
    p.add_argument("--mode", choices=list(MODES) + ["auto"], default="auto",
                   help="auto: 사이드체인이 있으면 heavy4.5, 없으면 ca8")
    p.add_argument("--compare", action="store_true",
                   help="계산하지 않고 --out CSV의 두 그룹을 나란히 요약한다")
    a = p.parse_args()

    if a.compare:
        compare(a.out)
        return
    if not (a.hotspots and a.group and a.pdbs):
        sys.exit("--hotspots, --group, PDB 가 모두 필요합니다 "
                 "(요약만 보려면 --compare)")

    want = []
    for tok in a.hotspots.split(","):
        tok = tok.strip()
        if tok:
            want.append((tok[0], int(tok[1:])))
    if not want:
        sys.exit("핫스팟이 비어 있습니다")

    files = []
    for g in a.pdbs:
        files.extend(sorted(glob.glob(g)))
    files = [f for f in files if os.path.isfile(f)]
    if not files:
        sys.exit("입력 PDB가 없습니다")

    # 모드 결정 — 첫 파일의 원자 구성을 보고 정한다
    if a.mode == "auto":
        _, has_sc = read_pdb(files[0], ca_only=False)
        mode = "heavy4.5" if has_sc else "ca8"
        print(f"[모드] auto -> {mode}  "
              f"({'사이드체인 있음' if has_sc else '백본만'})")
    else:
        mode = a.mode
    cfg = MODES[mode]
    print(f"[기준] {cfg['label']}, 핫스팟 {MIN_HIT}개 이상")

    rows, sources = [], defaultdict(int)
    for f in files:
        ch, _ = read_pdb(f, ca_only=cfg["ca_only"])
        if len(ch) != 2:
            print(f"[건너뜀] {os.path.basename(f)}: 사슬 {len(ch)}개", file=sys.stderr)
            continue
        # 타깃은 잔기 수가 많은 쪽(115) — binder는 70~100
        if a.target_chain:
            tc = a.target_chain
        else:
            tc = max(ch, key=lambda c: len(ch[c]))
        bc = [c for c in ch if c != tc][0]

        out2ref, src = target_numbering(f, tc, list(ch[tc]))
        sources[src] += 1

        binder_atoms = [xyz for atoms in ch[bc].values() for xyz in atoms]
        hit_out = contacts(binder_atoms, ch[tc], cfg["cutoff"])
        hit_ref = {out2ref.get(n, n) for n in hit_out}

        per = {f"{c}{n}": (n in hit_ref) for c, n in want}
        k = sum(per.values())
        rows.append(dict(
            name=os.path.splitext(os.path.basename(f))[0],
            group=a.group,
            binder_chain=bc, target_chain=tc,
            binder_len=len(ch[bc]),
            n_contact_total=len(hit_ref),
            n_hotspot_hit=k,
            success=int(k >= MIN_HIT),
            mode=mode,
            numbering_source=src,
            hotspots=";".join(f"{c}{n}" for c, n in want),
            **{f"hit{i+1}": int(per[f"{c}{n}"])
               for i, (c, n) in enumerate(want)},
        ))

    if not rows:
        sys.exit("계산된 설계가 없습니다")

    new = not os.path.isfile(a.out)
    with open(a.out, "a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        if new:
            w.writeheader()
        w.writerows(rows)

    n = len(rows)
    s = sum(r["success"] for r in rows)
    print(f"\n{a.group}: {n}건")
    print(f"  핫스팟 {MIN_HIT}개 이상 접촉: {s}/{n} = {s/n:.1%}   ← H2")
    print(f"  핫스팟 접촉 수 중앙값: "
          f"{sorted(r['n_hotspot_hit'] for r in rows)[n//2]} / {len(want)}")
    print(f"  계면 잔기 수 중앙값: "
          f"{sorted(r['n_contact_total'] for r in rows)[n//2]}")
    print("  핫스팟별 접촉률:")
    # 행의 키는 hit1..hit4 이고 want 순서와 같습니다. 잔기 이름으로 키를
    # 만들면 CSV 헤더가 먼저 쓰인 그룹의 것이어서 열이 어긋납니다.
    for i, (c, num) in enumerate(want, 1):
        print(f"    {c}{num}: {sum(r[f'hit{i}'] for r in rows)}/{n}")
    print(f"  번호 매핑 출처: {dict(sources)}")
    print(f"  -> {a.out}")
    if any(k != "trb" for k in sources):
        print("\n  [주의] trb 매핑을 쓰지 못한 설계가 있습니다. 번호 대응이")
        print("  맞는지 PDB를 직접 확인하고, 결과 보고 시 이 사실을 적으십시오.")


if __name__ == "__main__":
    main()
