# 자동 생성 — scripts/07_render.py. 직접 수정하지 마십시오.
# .pml 에는 레포 기준 상대 경로만 씁니다. 절대 경로를 박으면 커밋된
# 파일이 다른 사람의 기계에서 못 돕니다. **레포 루트에서 실행하십시오**
# (scripts/07_render.py 로 돌리면 자동으로 그렇게 됩니다).
reinitialize
set ray_opaque_background, 0
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set field_of_view, 20.0
# auto_zoom 을 끕니다 — 켜 두면 load 마다 카메라가 리셋되어 뒤에서 지정한
# set_view 가 중간 load 에 덮입니다 (레포 1 에서 겪은 함정).
set auto_zoom, off
set ray_interior_color, grey80
set two_sided_lighting, on
bg_color white

# Fig 1 — 4ZQK. PD-L1(사슬 A) 표면 위의 두 핫스팟 묶음.
#
# 답하는 질문: 대조군 패치가 **정말 다른 자리**인가.
# 두 패치 중심 사이 거리 25.8 A. 반대 면이라 한 장에 담기지 않아
# 앞/뒤 두 장으로 냅니다. 카메라는 계산된 set_view 이고 orient 가 아닙니다.
load inputs/PDL1_PD1_4ZQK.pdb, nat

hide everything
show surface, nat and chain A
color grey80, nat and chain A
color red, nat and chain A and resi 56+115+121+123
color marine, nat and chain A and resi 32+106+109+132

# native PD-1 — 실제 에피토프의 위치 표시. 얇게 둡니다.
show cartoon, nat and chain B
color palecyan, nat and chain B
set cartoon_transparency, 0.55, nat and chain B

# 대조 면 — 같은 분자를 반대쪽에서. PD-1 은 뒤로 갑니다
# 카메라 검증 (화면 좌표, 0 = 중앙 / ±1 = 가장자리). 눈으로 보지 않고
# 이 수치로 판정합니다 — 패치 중심은 반드시 (0.00, 0.00) 이어야 합니다.
#   패치 중심: 가로 -0.00  세로 -0.00
#   A32: 가로 -0.02  세로 +0.11
#   A106: 가로 -0.03  세로 -0.10
#   A109: 가로 +0.02  세로 -0.16
#   A132: 가로 +0.09  세로 +0.05
set_view (\
      0.5705,    -0.8213,     0.0000,\
     -0.5188,    -0.3603,     0.7753,\
     -0.6367,    -0.4423,    -0.6316,\
      0.0000,     0.0000,  -194.0825,\
     -9.8016,    45.3519,   114.2501,\
    166.9224,   220.9483,     0.0000 )
png figures/fig1b_control_face.png, width=1800, height=1350, dpi=300, ray=1
