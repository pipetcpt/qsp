# 혼합결합조직병 (MCTD) (Mixed Connective Tissue Disease, MCTD) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 자가면역·류마티스

[![MCTD QSP Model](mctd_qsp_model.png)](mctd_qsp_model.svg)

## 개요 (Overview)
혼합결합조직병(MCTD)은 전신 홍반 루푸스(SLE)·전신경화증(SSc)·다발성 근염(PM)의 임상 양상이 중복되어 나타나며, 항U1-RNP 항체 양성을 필수 혈청 표지자로 하는 자가면역 결합조직 질환입니다. 유병률은 10만 명당 약 3~7명으로 드물며, 레이노 현상이 거의 보편적으로(95%) 나타납니다. 핵심 병인은 항U1-RNP 항체 및 관련 T세포 면역 이상으로 촉발되는 다장기 염증이며, 폐동맥 고혈압(PAH)·간질성 폐 질환(ILD)·폐섬유화가 장기 예후를 결정하는 주요 합병증입니다. 치료는 스테로이드·하이드록시클로로퀸·MMF 기반으로 하며, PAH 합병 시 엔도텔린 수용체 길항제(보센탄 등)·PDE5 억제제를 추가합니다.

## 핵심 병태생리 경로 (Key Pathways)
| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| 항U1-RNP 항체 생성 | U1 snRNP에 대한 T세포·B세포 이중 활성화, IFN-α 신호 | 혈청 항U1-RNP 고역가 |
| T세포 이상 | Th1(IFN-γ)·Th17(IL-17) 증가, Treg 감소 | 다장기 염증 |
| 혈관 병변 | ET-1 증가·eNOS 감소·내피세포 손상 → 레이노·PAH | 폐동맥 혈관 리모델링 |
| 폐 섬유화 | TGF-β·IL-6에 의한 근섬유모세포 활성화 | FVC·DLCO 감소 |
| 근염 | IFN-γ 매개 근섬유 손상, CK 상승 | 근력 저하·MMT 감소 |
| 관절염 | TNF-α·IL-6·IL-17에 의한 활막 염증 | 관절통·관절염 |
| 보체 소모 | C3·C4 저하, 면역복합체 매개 보체 활성화 | SLE 양상 중복 |

## 주요 약물 표적 (Drug Targets)
- **하이드록시클로로퀸(HCQ)**: TLR7/9 억제 → IFN-α 생성 억제, 장기 기저 치료 (모든 환자)
- **프레드니솔론**: TNF-α·IL-6·IFN-γ 억제, 근염·관절염·장막염 조절
- **마이코페놀레이트(MMF/MPA)**: B세포·T세포 증식 억제 → 항체 생성 감소, ILD·신염
- **리툭시맙**: B세포 고갈, 중증 또는 난치성 경우
- **보센탄** (엔도텔린 수용체 길항제 ETA/ETB): ET-1 신호 차단 → 폐혈관 저항 감소 (PAH 병합 시)
- **타크로리무스**: 근염·신염 동반 시 칼시뉴린 억제로 T세포 조절

## 모델 구성 파일 (Model Files)
| 파일 | 설명 |
|------|------|
| [mctd_qsp_model.dot](mctd_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 182 노드 / 15 클러스터) |
| [mctd_qsp_model.svg](mctd_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [mctd_qsp_model.png](mctd_qsp_model.png) | PNG 이미지 (150 dpi) |
| [mctd_mrgsolve_model.R](mctd_mrgsolve_model.R) | mrgsolve ODE 모델 (약 38 구획 / 5개 치료 시나리오) |
| [mctd_shiny_app.R](mctd_shiny_app.R) | Shiny 대시보드 |
| [mctd_references.md](mctd_references.md) | 참고문헌 (약 35편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)
- **구획 구조**: 약물 PK 구획(HCQ 3구획, MPA/MMF 3구획 + 장간 순환, 프레드니솔론 2구획, 리툭시맙 2구획, 보센탄 2구획) + 질환 PD 구획(Th1, Th17, Treg, 나이브 B세포, GC B세포, 형질세포, 항U1-RNP 항체, TNF-α, IL-6, IL-17, IFN-γ, TGF-β, IFN-α, ET-1, 폐혈관 저항, 교원질/섬유화, FVC, DLCO, CK, MMT, 관절 종창 지수, C3, C4)
- **주요 치료 시나리오**: ① 자연경과, ② HCQ + 저용량 스테로이드, ③ HCQ + MMF + 스테로이드, ④ 리툭시맙 + MMF (중증), ⑤ 보센탄 추가 (PAH 합병 시 HCQ + 보센탄)
- **보정/근거**: EULAR/ERA-EDTA 결합조직병 동반 PAH 지침, SLE/SSc 주요 코호트(INPULSIS-MCTD 유사 연구), PARIS 레지스트리 데이터 정성적 참조

## Shiny 대시보드 (Dashboard)
환자 프로파일 입력(항U1-RNP 역가·기저 FVC·DLCO·PAP·근력·관절 종창 수), 약물 PK 농도 추이(HCQ·MMF·리툭시맙·보센탄), 주요 PD 바이오마커(사이토카인·항체·폐기능), 임상 엔드포인트(FVC/DLCO 변화·6분 보행 거리·MYOACT), 치료 시나리오 비교, 장기 폐 기능 예측 탭으로 구성됩니다.

## 실행 방법 (Usage)
```r
library(mrgsolve)
mod <- mread("mctd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("mctd_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg mctd_qsp_model.dot -o mctd_qsp_model.svg
```

## 참고문헌 (References)
자세한 인용은 [mctd_references.md](mctd_references.md) 참조 (약 35편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
