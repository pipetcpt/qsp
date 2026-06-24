# 중증 근무력증 (Myasthenia Gravis, MG) — QSP 모델

> **QSP 질환 모델 라이브러리** · Claude Code Routine으로 자동 생성된 정량적 시스템 약리학(QSP) 모델입니다.
> 상위 라이브러리 → [../README.md](../README.md) · 분류: 신경

[![MG QSP Model](mg_qsp_model.png)](mg_qsp_model.svg)

## 개요 (Overview)

중증 근무력증(MG)은 신경근접합부(NMJ)의 자가면역질환으로, 항아세틸콜린수용체(AChR) 항체(약 85%) 또는 항MuSK·항LRP4 항체에 의해 매개된다. 전 세계 유병률은 약 10만 명당 20명이며 여성과 젊은 연령대(20–30대)에서 초기 발병이 많고, 60세 이상 남성에서 두 번째 피크가 나타난다. 항체는 AChR을 직접 차단하거나 보체-매개 막공격복합체(MAC)를 형성하여 NMJ를 손상시키고, 이로 인해 변동성 근력 약화 및 피로가 특징적으로 발생한다. 흉선(흉선종 또는 흉선 과증식)에서 자가반응 T여포보조세포(TFH)가 GC B세포를 활성화하여 항체를 생산한다.

## 핵심 병태생리 경로 (Key Pathways)

| 경로 | 핵심 분자/기전 | 임상 결과 |
|------|----------------|-----------|
| TFH-GC B세포 축 | 흉선 TFH → GC B세포 → 단기·장기 형질세포 | 항AChR IgG 지속 생산 |
| AChR 항체 직접 차단 | IgG1/IgG3 → AChR 결합 및 cross-linking | 시냅스 후막 AChR 밀도 감소 |
| 보체 활성화 | C3b 옵소닌화 → MAC(C5b-9) 형성 | NMJ 구조 파괴, 근섬유 손상 |
| FcRn 매개 IgG 재활용 | FcRn이 IgG 분해에서 구출 → 항체 반감기 연장 | 지속성 고항체역가 |
| 콜린에스터분해효소 억제 | ACh 축적 → 더 많은 AChR 활성화 | 근력 일시 개선 |
| 흉선 병변 | 흉선종(10–15%) 또는 흉선 과증식 | 면역조절 이상 악화 |
| NMJ 안전계수 저하 | MEPP 감소 → 신경-근전달 실패 | QMG 점수 증가, 호흡 근력 저하 |

## 주요 약물 표적 (Drug Targets)

- **피리도스티그민 (Pyridostigmine)**: 아세틸콜린에스터분해효소(AChE) 억제 → ACh 축적, 증상 완화
- **프레드니솔론 (Prednisolone)**: 전반적 면역억제; TFH·B세포 활성 억제, Treg 확장
- **아자티오프린 (Azathioprine)**: 6-MP·6-TGN 전환을 통한 퓨린 합성 차단; 장기 유지요법
- **에쿨리주맙 (Eculizumab)**: C5 보체 억제 → MAC 형성 차단 (REGAIN 임상시험)
- **에프가르티기모드 (Efgartigimod)**: FcRn 차단 → IgG(항AChR 항체 포함) 혈중 농도 신속 감소
- **리툭시맙 (Rituximab)**: 항CD20 B세포 고갈; 항MuSK MG 및 난치성 MG에서 효과적

## 모델 구성 파일 (Model Files)

| 파일 | 설명 |
|------|------|
| [mg_qsp_model.dot](mg_qsp_model.dot) | Graphviz 기계론적 지도 소스 (약 202 노드 / 17 클러스터) |
| [mg_qsp_model.svg](mg_qsp_model.svg) | SVG 벡터 이미지 (확대 가능) |
| [mg_qsp_model.png](mg_qsp_model.png) | PNG 이미지 (150 dpi) |
| [mg_mrgsolve_model.R](mg_mrgsolve_model.R) | mrgsolve ODE 모델 (약 20 구획 / 7개 치료 시나리오) |
| [mg_shiny_app.R](mg_shiny_app.R) | Shiny 대시보드 |
| [mg_references.md](mg_references.md) | 참고문헌 (약 41편, PubMed 링크) |

## mrgsolve 모델 (ODE Model)

- **구획 구조**: 피리도스티그민 경구 2구획 PK, 프레드니솔론 2구획, 아자티오프린(GUT→6-MP→6-TGN) 3구획, 에쿨리주맙 2구획+C5 결합, 에프가르티기모드 1구획+FcRn, 흉선 TFH·GC B세포·SLPC·LLPC 세포 구획, 항AChR IgG, AChR 밀도, 보체 활성화, NMJ 안전계수 및 QMG 점수
- **주요 치료 시나리오**: ① 무치료 MG ② 피리도스티그민 단독 ③ 피리도스티그민+프레드니솔론 ④ 피리도스티그민+프레드니솔론+아자티오프린 ⑤ 에쿨리주맙+피리도스티그민 ⑥ 에프가르티기모드+피리도스티그민 ⑦ 리툭시맙+피리도스티그민
- **보정/근거**: REGAIN(에쿨리주맙), ADAPT(에프가르티기모드), MGTX(흉선절제술) 임상시험 데이터 참조

## Shiny 대시보드 (Dashboard)

6개 탭으로 구성: ① **환자 프로파일** (항체 유형·기저 QMG·흉선 상태 설정), ② **PK** (혈장 약물 농도, FcRn 점유율, 에쿨리주맙 C5 억제율), ③ **PD 주요지표** (항AChR IgG·AChR 밀도·보체 활성화 추이), ④ **임상 엔드포인트** (QMG, MGC, NMJ 안전계수), ⑤ **시나리오 비교** (7개 치료 전략 직접 비교), ⑥ **바이오마커** (IgG 감소율, SLPC/LLPC 변화).

## 실행 방법 (Usage)

```r
library(mrgsolve)
mod <- mread("mg_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny 대시보드:
shiny::runApp("mg_shiny_app.R")
```
```bash
# 기계론적 지도 렌더링
dot -Tsvg mg_qsp_model.dot -o mg_qsp_model.svg
```

## 참고문헌 (References)

자세한 인용은 [mg_references.md](mg_references.md) 참조 (약 41편).

---
*본 모델은 교육·연구 목적의 정성적·반정량적 QSP 모델이며, 실제 임상 의사결정에 직접 사용할 수 없습니다.*
