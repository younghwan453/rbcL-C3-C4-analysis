# 11주차 · 생물학적 배경과 연구 설계

## 목표와 수행 내용

광합성 경로와 Rubisco의 기능을 이해하고, 문헌에 따라 네 그룹에서 8·7·8·7종을 선정했습니다. 11주차 초기 목록과 속·과 비교 근거는 [후보종 설명](species_selection.md)에, 현재의 검증 상태와 종별 출처는 [메타데이터 CSV](../data/processed/species_metadata.csv)에 정리했습니다.

[후보 검토 노트북](analysis.ipynb)에서 그룹별 집계, 시험 수집 대상, 속·과별 대비를 재현할 수 있습니다.

| 분석 그룹 | 후보 수 |
| --- | ---: |
| C3 외떡잎식물 | 8 |
| C3 진정쌍떡잎식물 | 7 |
| C4 외떡잎식물 | 8 |
| C4 진정쌍떡잎식물 | 7 |

## 사용한 프로그램과 기록

PowerShell, Git, GitHub CLI, Conda와 Python을 사용해 실행 환경과 변경 사항을 확인합니다. [연구 계획](research_plan.md), [학습 자료](biology_background.md), [학생 학습 기록](learning_log.md), [AI 협업 기록](ai_usage_log.md), [손 계산 예제](hand_calculation.md), [테스트 기록](test_results.md), [참고문헌](references.md)을 함께 검토할 수 있습니다.

## 생물학적 해석, 한계와 다음 주차

같은 속의 Flaveria, Atriplex, Cyperus와 같은 과의 Cleomaceae에서 C3/C4 대비를 구성했습니다. 이러한 대비는 분류 범위가 넓은 전체 비교를 보완합니다. 가까운 종의 서열은 공통조상 때문에 서로 비슷할 수 있으며, 표본 할당에서 C3/C4와 분류군의 비율을 맞췄으므로 모집단 비율을 추정할 수 없습니다.

11주차 결과는 후보 목록과 분석 설계입니다. 생물학적 서열 결과와 그래프의 생성 시점은 12주차부터입니다. 다음 주에는 그룹별로 Oryza sativa, Flaveria pringlei, Zea mays, Flaveria bidentis를 먼저 시험 수집하고, 네 레코드의 검증이 끝난 뒤 30종에 적용합니다.
